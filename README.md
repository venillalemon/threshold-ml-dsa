# threshold-ml-dsa

n 方恶意安全的 **ML-DSA (FIPS 204) 门限签名** MPC 实现:私钥以秘密份额形式
分散在 n 个参与方手中,任何一方都拿不到完整私钥,各方合作产出一个**标准的
ML-DSA 签名**(验证方无需知道签名是 MPC 产的)。

协议来自论文的 Π_MLDSA(KeyGen + Sign),实现为单次执行(T=1)变体:
拒绝采样失败时返回 (c, ⊥, ⊥),重跑一次 Sign 即可(单次通过率约 20%~35%,
取决于参数集)。

三个计算域协同工作(**混合安全后端**,见下):

- **算术域 F_q**(q = 8380417):BDOZ 成对认证加性份额,负责线性代数
  (A·y + e 这类矩阵乘,本地免费),开启一律"收到即验、单 flight";
- **离线布尔域(C_pre)**:**认证 GMW / TinyOT**(真实恶意 COT 三元组),
  跑挑战前的电路——A2B 恢复 w、Decompose、边界 producer;
- **在线布尔域(C_post)**:**WRK17b 四行认证混淆**,在**离线阶段**就把电路
  混淆好、四行表和固定输入标签发给求值方 P1、并把**输出掩码也在离线开给 P1**,
  于是挑战后只剩**一次 δ_H 开启 + 一次标签投递 + 本地解码**(两个 flight)。
- 三域之间靠 **edaBits/daBits**(同一值的算术 + GMW 布尔一致份额对)和
  **A2B 转换**(mask-and-open)过桥。

## 混合安全后端(GMW 离线 / WRK 在线)

早期原型把所有布尔电路都放在单一 KRRW 认证混淆会话里,求值后需要广播
Λ_AND、掷币、commit-open 零检查等多轮交互,在线轮数远超论文目标。现在按论文
的 offline/online 切分重构:

| 组件 | 文件 | 作用 |
|---|---|---|
| `src/infra/backend.h` | `Backend<nP>` | 一个 `NetIOMP` 同时供 GMW(离线电路)、WRK(在线电路)与算术 BDOZ 开启;每方私有 pinned Δ |
| `src/circuit/wrk_phase2.h` | `wrk_offline` / `wrk_online` | C_post 离线混淆 + 离线开输出掩码;在线一次标签投递后本地解码 |
| `src/protocol/sign_2round.h` | `cpre_program` | 两轮:**一张** GMW 电路 C_pre(恢复 w、Decompose、producer);C_post 的 late 输入是**公开的** δ_H(`wrk_online`) |
| `src/protocol/sign_slot.h` | `cpre_program` | slot:**一张** GMW 电路 C_pre(恢复 w、Decompose、producer、one-hot 解码、压缩网络、空槽 MUX、ovf);C_post 的 late 输入 d^ 由 P1 **free-XOR** u 的标签得到(`wrk_online_routed`) |
| `third_party/emp-ag/` | `emp::gmw` / `emp::wrk` | vendored WRK 四行 + 认证 GMW 后端(依赖 emp-tool / emp-ot) |

**状态**:两种模式都在新后端上跑通并验证(2/3 方 ML-DSA-44):KeyGen/Decompose
自检 0 错,电路判定与明文一致,产出的 (c,z,h) 通过真实 ML-DSA 验证,TAMPER_C 在
checked open 处 abort,**在线都是 2 个 flight**。

slot 模式的 C_post 里,d^_j[ℓ] = ⊕_{i: δ_H,i[ℓ]=1} u_{j,i} 由 P1 对 u 的活标签做
free-XOR 得到,零个门、无 PublicBits(与 slot 论文 Fig. 1 一致)。与论文文字唯一
的出入:d^ 后面每根线接一个 `AND 1` 缓冲门,把依赖路由的掩码换成新鲜掩码——这
NS·ν = 576 个门(ML-DSA-44)的表只能在 flight 1 之后由混淆方本地造好、随 flight 2
一起发给 P1;其余(加法器、AND 树、r·bin23(ρ))全部离线混淆。代价是 online 门数
24,127 → 24,703、n=3 在线通信 0.27 → 0.38 MB,轮数仍是 2。

## 仓库结构

```
threshold-ml-dsa/
├── main.cpp              入口:建 session 和 FakeDealer,跑 KeyGen + Sign
├── src/
│   ├── infra/            基础设施(demo 作弊集中在这一层)
│   ├── circuit/          布尔电路组件
│   └── protocol/         协议——每个文件对应论文一张协议图
├── docs/                 Decompose 电路门数推导(decompose-proof.tex)
├── runN                  n 方本地启动脚本(make run 依赖它)
└── Makefile / CMakeLists.txt
```

### src/infra/ — 基础设施

| 文件 | 内容 |
|---|---|
| `ref.h` | ML-DSA 参数常量、F_q/R_q 运算、`ref_*` 明文参考实现(FIPS 204)、ExpandA/H/SampleInBall 的 **stub** |
| `bdoz.h` | **真实协议层**:成对认证(BDOZ)的 `FqShare`/`RingShare`、收到即验的 checked open、`matvec_negacyclic` |
| `dealer.h` | **所有 demo 作弊集中一处**:`FakeDealer` 从共享种子发 α、Δ,`deal_fq/ring/bits` 原语,edaBit 家族 `deal_y_edabit`/`deal_fq_edabit`/`deal_ring_edabit`、小范数 `deal_secret_poly`;`SharePair`/`RingEdabits` 类型 |
| `backend.h` | `Backend<nP>`:一个 `NetIOMP` 供 GMW/WRK/算术开启;Δ 默认私有随机,demo 传入 dealer 的 Δ |

### src/circuit/ — 布尔电路

| 文件 | 内容 |
|---|---|
| `circuit.h` | 通用小组件:fa/ha 加法器、行波加法/常数加、有符号比较、or_tree、pack/unpack |
| `decompose.h` | FIPS Decompose 电路(43 AND / 100 AND 每系数,推导见 `docs/`) |
| `a2b.h` | F_A2B:mask-and-open 把 `<x>_q` 转成电路内的 `<x>_2`;含 `sub_modq` 模减电路 |

### src/protocol/ — 协议层

| 文件 | 对应论文 | 内容 |
|---|---|---|
| `keygen.h` | Π_MLDSA.KeyGen | dealer 采 s/e 并发环份额,**明文**算 t = A·s+e,Power2Round,产出 pk 与 `KeyPair`(零网络) |
| `sign_2round.h` | Π_TwoRound(main.pdf) | 步骤 1 eDaBits(dealer)→ 步骤 2 w = Ay+e_w、开 c(BDOZ)→ 步骤 3 一张 C_pre(GMW)→ C_post 离线混淆(WRK,δ_H 公开 late 输入)→ 挑战、flight 1 开 δ_H、flight 2 投标签、本地解码、MakeHint |
| `sign_slot.h` | Π_Slot(slot pdf Fig. 1) | 同上结构;C_pre 多了解码/压缩/MUX,开 (w1, ovf);C_post 固定输入 M̂、ρ、u,d^ 由 free-XOR 路由,flight 2 携带 N_s·ν 个缓冲门的表 |

分层依赖自下而上:`infra → circuit → protocol → main.cpp`。

## 当前的 demo 作弊(与真实实现的差距)

协议逻辑忠实于论文;以下基础原语是占位实现,替换时协议层不需要改动:

1. **edaBits/daBits**(`infra/dealer.h`):dealer 直接在正确范围内采样(y、
   R2、R_H、以及小范数 s/e/e_w),一次发出该值的算术(F_q + 环)与 GMW/WRK
   布尔三种一致份额——**零通信、零 GMW、无拒绝采样**。真实现需在不公开值的
   前提下产出份额对(如 MP-SPDZ 的 edaBits 协议)。因此 KeyGen 本身不跑任何
   电路(GMW ANDs = 0);GMW 只用于 Sign 的实际电路(A2B、Decompose、
   producer)。
2. **SPDZ / GMW offline**(`infra/dealer.h` 的 `FakeDealer`):全量算术 α 与每方
   GMW/WRK 的 Δ 都由共享种子确定性导出,dealer 借此伪造一致的布尔认证份额。
   真实现中 α、Δ 从不被任何一方重构,认证份额由 MASCOT/Overdrive(算术)与
   真实恶意 COT(布尔)产出;协议层只接触 `deal_*` 和 `my_alpha`/`my_delta`,
   替换面已隔离;
3. **MAC 域太小**:MAC 在 F_q(2²³)上,伪造概率 2⁻²³;真实现需扩域或多重 MAC;
4. **SHAKE 全家**(`infra/ref.h`):ExpandA/H/SampleInBall 用 seed_seq+mt19937
   代替 SHAKE-128/256,均为公开本地计算,换真 Keccak 不改协议消息。

## 依赖与安装

- C++20 编译器、CMake ≥ 3.25、OpenSSL
- **emp-tool** 和 **emp-ot** `cmake --install` 到标准前缀(默认 /usr/local):
  - [emp-tool](https://github.com/emp-toolkit/emp-tool)(基础:NetIO、block、PRG、电路 IR)
  - [emp-ot](https://github.com/emp-toolkit/emp-ot)(OT 扩展:Ferret / SoftSpoken / IKNP)
- **WRK/GMW 后端已 vendored 在 `third_party/emp-ag/`**,由本仓库的
  CMake `add_subdirectory` 带出 `emp-ag::gmw`(自动依赖 emp-tool / emp-ot),
  无需单独安装。

emp-tool / emp-ot 安装方式相同:

```bash
cmake -S . -B build && cmake --build build -j && sudo cmake --install build
```

## 跑

```bash
make n=3              # 3 方恶意,编译并本机跑一次 KeyGen + Sign
make n=5              # 5 方(n= 任意值,nP 是模板参数,现编对应二进制)
make n=3 p=44         # 换 ML-DSA-44 参数集(默认 65;还有 p=87)
make n=3 TEST=1       # 测试构建:开出秘密值做范围/电路-明文逐位比对/签名验证
make n=3 TEST=1 TAMPER_C=1   # 对抗测试:party 1 篡改一个 share,MACCheck 应 abort
make clean
```

默认 slot 模式(`SLOT=1`),`SLOT=0` 切两轮基线。
端口每次随机(避开上一轮 TIME_WAIT);要固定就 `PORT=16400`。

后端自身还有独立联网自检(2/3/5 方的 GMW COT+open、`gmw.evaluate`、
`wrk_offline`/`wrk_online`),见 `tests/`:

```bash
cmake -S . -B build -DMLDSA_TESTS=ON -DMLDSA_DRIVER=OFF && cmake --build build -j
# 每方一进程,共享 EMP_PORT,例如 3 方:
for p in 1 2; do EMP_PORT=17720 ./build/wrk_bridge_test3 $p & done; EMP_PORT=17720 ./build/wrk_bridge_test3 3
```

成功输出(绿色为完整签名;红色 `(c, bot, bot)` 是 T=1 下正常的拒绝采样,重跑即可):

```
  pk: rho[0]=9e3779b9, t1=1536 coefficients | tr[0]=...
  signature: (c, z, h) — ||z||_inf=130731 (< 130994), HW(h)=67 (<= 80)
  setup     42 ms | comm    0.5 MB
  KeyGen    15 ms | comm    0.0 MB | GMW ANDs 0
  Sign    1456 ms | comm   77.8 MB | GMW ANDs 398057 | G2(WRK) ANDs 78079
    - offline (no msg)   1425 ms | comm   77.3 MB
    - online  (after msg)  30 ms | comm    0.5 MB | flights 2
STATS np=2 param=ML-DSA-65 slot=0 ...
  KeyGen + Sign complete
```

时间/通信量按 **离线 / 在线** 拆分:重活(GMW 三元组、C_post 混淆)全在离线,
在线只有 2 个 flight。`STATS` 行是机器可读版,供基准脚本 grep。

TEST=1 会额外打印每步自检:KeyGen 的 t = A·s+e 明文比对与 Power2Round 恒等式、
PrepSign 的 Decompose 对 FIPS oracle、Sign 的 r₀/z 电路对明文逐位比对,
以及产出完整签名时跑一遍真验证(UseHint → w1′ → 重算挑战)。

### runN 与手动多终端

`make n=3` 底层调 `./runN 3 ./build/main3`:party 1..n−1 错峰后台启动、
party n 前台,并聚合所有方的退出码(任何一方 MACCheck abort 都会让 make 报错)。

调试时可手动分终端跑(`argv[1]` 是 party 编号,`EMP_PORT` 各方一致):

```bash
EMP_PORT=16400 ./build/main3 1    # 自检输出在 party 1
EMP_PORT=16400 ./build/main3 2
EMP_PORT=16400 ./build/main3 3
```
