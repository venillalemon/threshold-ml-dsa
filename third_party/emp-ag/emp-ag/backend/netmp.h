#ifndef EMP_AG_BACKEND_NETMP_H__
#define EMP_AG_BACKEND_NETMP_H__
#include "emp-tool/runtime/io/net_io_channel.h"
#include <array>
#include <atomic>
#include <arpa/inet.h>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <netinet/in.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace emp {

template <int nP>
inline constexpr bool valid_netmp_party(int party) {
  return party >= 1 && party <= nP;
}

template <int nP>
inline constexpr bool valid_netmp_base_port(int port) {
  return port >= 1 && (int64_t)port + nP <= 65535;
}

template <int nP>
inline constexpr bool valid_netmp_mesh_count(int n_meshes) {
  return n_meshes >= 1 && n_meshes <= (1 << 16) &&
         (int64_t)n_meshes * nP <= INT_MAX;
}

// ---- mesh connection setup core -------------------------------------------
// ONE listening socket per party (`port + party`) serves EVERY mesh: each
// dialed connection announces {party, mesh_tag} (4 bytes, written before any
// protocol byte and before the FS transcript hash is enabled, so the handshake
// never enters a digest), and the accept loop routes by that tag. Arrival
// ORDER therefore never matters — a fast peer may already be dialing mesh 1
// while we still wait on mesh-0 accepts; the tag, not the accept queue,
// decides where a connection lands. This is what lets k meshes share one
// port: k separate NetIOMP ctors on offset ports would multiply the port
// range (and the old `port + 100` convention collided outright at nP >= 100),
// while sequential same-port reuse races the kernel accept queue.
//
// Fills dialed[m][peer] (connections THIS party initiated) and
// accepted[m][peer] (connections the peer initiated), 1-based by peer,
// [party] slots null. Out-of-range ids, bad tags, and duplicate (tag, peer)
// pairs abort before they can index or overwrite anything.
template <int nP>
inline void netmp_setup(int party, int port, const char *const *ip, int n_meshes,
                        std::vector<std::array<NetIO *, nP + 1>> &dialed,
                        std::vector<std::array<NetIO *, nP + 1>> &accepted) {
  if (!valid_netmp_party<nP>(party)) {
    fprintf(stderr, "netiomp: party must be in [1, nP]\n"); exit(1);
  }
  if (!valid_netmp_base_port<nP>(port)) {
    fprintf(stderr, "netiomp: base port must keep port+nP in [1, 65535]\n"); exit(1);
  }
  if (!valid_netmp_mesh_count<nP>(n_meshes)) {
    fprintf(stderr, "netiomp: mesh count must be in [1, 65536] and fit setup arithmetic\n");
    exit(1);
  }
  dialed.assign(n_meshes, {});
  accepted.assign(n_meshes, {});

  int listen_fd = socket(AF_INET, SOCK_STREAM, 0);
  if (listen_fd < 0) { perror("netiomp socket"); exit(1); }
  int reuse = 1;
  setsockopt(listen_fd, SOL_SOCKET, SO_REUSEADDR, &reuse, sizeof(reuse));
  sockaddr_in serv{};
  serv.sin_family = AF_INET;
  serv.sin_addr.s_addr = htonl(INADDR_ANY);
  serv.sin_port = htons(port + party);
  // Bounded retry on EADDRINUSE, mirroring the connect loop below: back-to-back
  // localhost runs can land on a port whose previous socket the kernel is still
  // tearing down (SO_REUSEADDR clears a TIME_WAIT peer but not a socket in the
  // last moments of close). Wait that out rather than aborting a whole party on
  // the first collision; any other bind error is a real misconfiguration and
  // fails immediately.
  for (int tries = 0; ; ++tries) {
    if (bind(listen_fd, (sockaddr *)&serv, sizeof(serv)) == 0) break;
    if (errno != EADDRINUSE || tries >= 5000) { perror("netiomp bind"); exit(1); }
    usleep(1000);   // ~5 s ceiling; a released port frees in well under that
  }
  if (listen(listen_fd, n_meshes * nP) < 0) { perror("netiomp listen"); exit(1); }

  auto connect_to = [&](int peer, uint16_t tag) -> int {
    int peer_port = port + peer;
    // Validate the peer address once: inet_addr returns INADDR_NONE on a malformed
    // string, which would otherwise be sent as the (broadcast) 255.255.255.255.
    const char *host = ip ? ip[peer] : "127.0.0.1";
    in_addr_t addr = inet_addr(host);
    if (addr == INADDR_NONE && strcmp(host, "255.255.255.255") != 0) {
      fprintf(stderr, "netiomp: invalid peer %d address '%s'\n", peer, host); exit(1);
    }
    // Bounded retry: a peer may still be starting up, but never spin forever — an
    // unbounded loop hangs the whole mesh setup silently if the peer never arrives.
    for (int tries = 0; tries < 120000; ++tries) {   // ~120 s at 1 ms/try
      int sock = socket(AF_INET, SOCK_STREAM, 0);
      if (sock < 0) { perror("netiomp socket"); exit(1); }
      sockaddr_in dest{};
      dest.sin_family = AF_INET;
      dest.sin_addr.s_addr = addr;
      dest.sin_port = htons(peer_port);
      if (connect(sock, (sockaddr *)&dest, sizeof(dest)) == 0) {
        uint16_t hs[2] = {(uint16_t)party, tag};
        if (write(sock, hs, sizeof(hs)) != (ssize_t)sizeof(hs)) {
          fprintf(stderr, "netiomp handshake send failed\n"); exit(1);
        }
        return sock;
      }
      close(sock);
      usleep(1000);
    }
    fprintf(stderr, "netiomp: could not connect to peer %d after ~120 s\n", peer); exit(1);
  };

  // One connector thread per peer, dialing every mesh in tag order; each
  // thread writes only its own peer's dialed[*][peer] slots (race-free).
  std::vector<std::thread> connectors;
  connectors.reserve(nP - 1);
  for (int p = 1; p <= nP; ++p) if (p != party) {
    connectors.emplace_back([&, p] {
      for (int m = 0; m < n_meshes; ++m) {
        int sock = connect_to(p, (uint16_t)m);
        dialed[m][p] = new NetIO(sock, true);
      }
    });
  }

  // Main thread accepts the matching n_meshes * (nP - 1) connections, in
  // whatever order they arrive, and routes each by its handshake.
  for (int k = 0; k < n_meshes * (nP - 1); ++k) {
    sockaddr_in cli; socklen_t clilen = sizeof(cli);
    int fd = accept(listen_fd, (sockaddr *)&cli, &clilen);
    if (fd < 0) { perror("netiomp accept"); exit(1); }
    // Bound the handshake read so a peer that connects but never sends its id
    // can't stall mesh setup forever; reset to blocking once the 4-byte id
    // arrives (the fd then carries the whole protocol with no timeout).
    timeval tv{30, 0};
    setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &tv, sizeof(tv));
    uint16_t hs[2];
    if (read(fd, hs, sizeof(hs)) != (ssize_t)sizeof(hs)) {
      fprintf(stderr, "netiomp handshake recv failed or timed out\n"); exit(1);
    }
    timeval none{0, 0};
    setsockopt(fd, SOL_SOCKET, SO_RCVTIMEO, &none, sizeof(none));
    const int peer = (int)hs[0], tag = (int)hs[1];
    // Reject a bogus handshake before it indexes dialed/accepted (a valid peer is
    // a party in 1..nP other than ourselves; a valid tag names a mesh) — otherwise
    // it is an OOB write. A duplicate (tag, peer) would silently orphan a channel.
    if (peer < 1 || peer > nP || peer == party) {
      fprintf(stderr, "netiomp handshake: invalid peer id %d\n", peer); exit(1);
    }
    if (tag < 0 || tag >= n_meshes) {
      fprintf(stderr, "netiomp handshake: invalid mesh tag %d\n", tag); exit(1);
    }
    if (accepted[tag][peer] != nullptr) {
      fprintf(stderr, "netiomp handshake: duplicate connection from peer %d mesh %d\n",
              peer, tag); exit(1);
    }
    accepted[tag][peer] = new NetIO(fd, true);
  }
  for (auto &t : connectors) t.join();

  close(listen_fd);
}

template <int nP> class NetIOMP {
public:
  NetIO *ios[nP + 1];
  NetIO *ios2[nP + 1];
  int party;

  // Symmetric pairwise mesh. For each peer p != party there are two TCP
  // connections so each direction has its own duplex channel. Convention:
  // for pair (i, j) with i < j, the smaller party DIALS the `ios` slot's
  // connection; the larger party dials the `ios2` slot's. Setup runs through
  // netmp_setup (one listening port per party; see above). `ip` supplies a
  // per-party address table (1-based; ip[p] for peer p); nullptr means
  // loopback for all. For several meshes over the SAME port, use NetIOMPGroup.
  NetIOMP(int party, int port, const char *const *ip = nullptr) {
    std::vector<std::array<NetIO *, nP + 1>> dialed, accepted;
    netmp_setup<nP>(party, port, ip, /*n_meshes=*/1, dialed, accepted);
    adopt_(party, dialed[0], accepted[0]);
  }

  // Adopt pre-connected channels (used by NetIOMPGroup and the ctor above):
  // routes this party's dialed / accepted connections for one mesh into the
  // ios/ios2 slots per the smaller-party-dials-`ios` convention. Takes
  // ownership (the dtor deletes every channel).
  NetIOMP(int party, std::array<NetIO *, nP + 1> &dialed,
          std::array<NetIO *, nP + 1> &accepted) {
    adopt_(party, dialed, accepted);
  }

  int64_t count() {
    int64_t res = 0;
    for (int i = 1; i <= nP; ++i)
      if (i != party) {
        res += ios[i]->send_counter + ios[i]->recv_counter;
        res += ios2[i]->send_counter + ios2[i]->recv_counter;
      }
    return res;
  }

  ~NetIOMP() {
    for (int i = 1; i <= nP; ++i)
      if (i != party) {
        delete ios[i];
        delete ios2[i];
      }
  }

  // Owns the NetIO* channels it deletes above; non-copyable / non-movable so the
  // raw owning pointers are never double-deleted.
  NetIOMP(const NetIOMP &) = delete;
  NetIOMP &operator=(const NetIOMP &) = delete;
  NetIOMP(NetIOMP &&) = delete;
  NetIOMP &operator=(NetIOMP &&) = delete;

  void send_data(int dst, const void *data, size_t len) {
    if (dst != 0 and dst != party) {
      if (party < dst)
        ios[dst]->send_data(data, len);
      else
        ios2[dst]->send_data(data, len);
    }
  }
  void send_bool(int dst, const bool *data, int64_t len) {
    if (dst != 0 and dst != party) {
      if (party < dst)
        ios[dst]->send_bool(data, len);
      else
        ios2[dst]->send_bool(data, len);
    }
  }
  void recv_data(int src, void *data, size_t len) {
    if (src != 0 and src != party) {
      // No auto-flush of the send stream here: callers already flush
      // before any round-trip, and acquiring flockfile() on the send
      // stream from a recv thread can deadlock against a pool thread
      // that's blocked in write() holding that same stdio lock.
      if (src < party)
        ios[src]->recv_data(data, len);
      else
        ios2[src]->recv_data(data, len);
    }
  }
  void recv_bool(int src, bool *data, int64_t len) {
    if (src != 0 and src != party) {
      if (src < party)
        ios[src]->recv_bool(data, len);
      else
        ios2[src]->recv_bool(data, len);
    }
  }
  NetIO *&get(size_t idx, bool b = false) {
    if (b)
      return ios[idx];
    else
      return ios2[idx];
  }
  void flush(int idx = 0) {
    if (idx == 0) {
      for (int i = 1; i <= nP; ++i)
        if (i != party) {
          ios[i]->flush();
          ios2[i]->flush();
        }
    } else {
      if (party < idx)
        ios[idx]->flush();
      else
        ios2[idx]->flush();
    }
  }

private:
  void adopt_(int party_, std::array<NetIO *, nP + 1> &dialed,
              std::array<NetIO *, nP + 1> &accepted) {
    party = party_;
    for (int i = 0; i <= nP; ++i) { ios[i] = nullptr; ios2[i] = nullptr; }
    for (int p = 1; p <= nP; ++p) if (p != party) {
      if (party < p) {          // we are the smaller party of the pair
        ios[p]  = dialed[p];    //   we dialed the ios connection
        ios2[p] = accepted[p];  //   the peer dialed ios2
      } else {
        ios[p]  = accepted[p];  //   the (smaller) peer dialed ios
        ios2[p] = dialed[p];    //   we dialed ios2
      }
    }
  }
};

// Several INDEPENDENT NetIOMP meshes over ONE listening port per party: all
// k * (nP - 1) inbound connections arrive on the same port and are routed by
// the {party, mesh_tag} handshake (netmp_setup, above), so adding a mesh —
// e.g. a dedicated preprocessing mesh for the cross-chunk overlap — adds
// sockets but NO ports, no offset arithmetic, and no accept-order races.
// Every party must construct its group with the same n_meshes.
template <int nP>
class NetIOMPGroup {
public:
  NetIOMPGroup(int party, int port, int n_meshes, const char *const *ip = nullptr) {
    std::vector<std::array<NetIO *, nP + 1>> dialed, accepted;
    netmp_setup<nP>(party, port, ip, n_meshes, dialed, accepted);
    meshes_.reserve(n_meshes);
    for (int m = 0; m < n_meshes; ++m)
      meshes_.emplace_back(std::make_unique<NetIOMP<nP>>(party, dialed[m], accepted[m]));
  }
  NetIOMP<nP> &mesh(int m) { return *meshes_[(size_t)m]; }
  const NetIOMP<nP> &mesh(int m) const { return *meshes_[(size_t)m]; }
  int size() const { return (int)meshes_.size(); }

  NetIOMPGroup(const NetIOMPGroup &) = delete;
  NetIOMPGroup &operator=(const NetIOMPGroup &) = delete;

private:
  std::vector<std::unique_ptr<NetIOMP<nP>>> meshes_;
};

}  // namespace emp
#endif // EMP_AG_BACKEND_NETMP_H__
