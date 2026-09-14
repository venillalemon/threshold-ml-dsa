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
| `src/protocol/prepsign.h` | `recover_decompose_program` | A2B + Decompose 编译成一张 GMW 电路 |
| `src/protocol/sign_2round.h` | `producer_program` | 边界 producer 编译成一张 GMW 电路;C_post 走 WRK 桥 |
| `third_party/emp-ag/` | `emp::gmw` / `emp::wrk` | vendored WRK 四行 + 认证 GMW 后端(依赖 emp-tool / emp-ot) |

**状态**:两轮模式(`SLOT=0`,即论文 main.pdf 实现)已完整跑通并验证——
2/3 方 ML-DSA-44 下 KeyGen/Decompose 自检 0 错,电路判定与明文一致,产出的
(c,z,h) 通过真实 ML-DSA 验证;TAMPER_C 在 checked open 处 abort;在线 2 个
flight。**slot 模式(`SLOT=1`,sign.h)尚未迁移**到新后端(它的 online 路由
输入 Dh 是"秘密 late 输入",需要给 WRK 桥加一条在线 masked-open 安装路径)。

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
| `spdz.h` | SPDZ 认证份额 `AuthShare`、open、MACCheck(SPDZ-2 式)、`FakeDealer`(全量 α 的唯一居所) |
| `edabits.h` | `SharePair` 类型;unsigned/signed/y/Fq 四种 edaBits(**demo cheat**:开 r 后 dealer 拆分) |

### src/circuit/ — 布尔电路

| 文件 | 内容 |
|---|---|
| `circuit.h` | 通用小组件:fa/ha 加法器、行波加法/常数加、有符号比较、or_tree、pack/unpack |
| `decompose.h` | FIPS Decompose 电路(43 AND / 100 AND 每系数,推导见 `docs/`) |
| `a2b.h` | F_A2B:mask-and-open 把 `<x>_q` 转成电路内的 `<x>_2`;含 `sub_modq` 模减电路 |

### src/protocol/ — 协议层

| 文件 | 对应论文 | 内容 |
|---|---|---|
| `rand.h` | Π_smallnormpoly | daBit 拒绝采样生成 [-η,η] 均匀小范数多项式份额 |
| `keygen.h` | Π_MLDSA.KeyGen | s/e 采样、开 t + MACCheck、Power2Round,产出 pk 与 `KeyPair` |
| `prepsign.h` | Π_PrepSign | 离线预处理:采 y/e_w、w = Ay+e_w、A2B、Decompose,返回 (⟨w₀⟩₂, w₁, ⟨y⟩₂) |
| `sign.h` | Π_MLDSA.Sign | 在线签名(T=1):挑战 c、r₀/z 拒绝电路(只公开判定位)、开 z、MakeHint |

分层依赖自下而上:`infra → circuit → protocol → main.cpp`。

## 当前的 demo 作弊(与真实实现的差距)

协议逻辑忠实于论文;以下基础原语是占位实现,替换时协议层不需要改动:

1. **edaBits/daBits**(`infra/edabits.h`):布尔份额抽好后把 r 公开、再由
   dealer 拆出一致的算术份额。真实现需在不公开 r 的前提下产出份额对
   (如 MP-SPDZ 的 edaBits 协议);
2. **SPDZ offline**(`infra/spdz.h` 的 `FakeDealer`):全量 MAC 钥匙 α 由共享
   种子导出。真 SPDZ 中 α 从不被任何一方重构,认证份额由 MASCOT(OT)或
   Overdrive(HE)产出;协议层只接触 `my_alpha` 与 `deal()`,替换面已隔离;
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

默认走两轮模式(`SLOT=0`);slot 模式尚未迁移(见上「混合安全后端」)。
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
  KeyGen   173 ms | comm    3.6 MB | GMW ANDs 17424
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
