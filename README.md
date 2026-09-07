# threshold-ml-dsa

n 方恶意安全的 **ML-DSA (FIPS 204) 门限签名** MPC 实现:私钥以秘密份额形式
分散在 n 个参与方手中,任何一方都拿不到完整私钥,各方合作产出一个**标准的
ML-DSA 签名**(验证方无需知道签名是 MPC 产的)。

协议来自论文的 Π_MLDSA(KeyGen + Sign),实现为单次执行(T=1)变体:
拒绝采样失败时返回 (c, ⊥, ⊥),重跑一次 Sign 即可(单次通过率约 20%~35%,
取决于参数集)。

两个计算域协同工作:

- **算术域 F_q**(q = 8380417):SPDZ 式认证加性份额,负责线性代数
  (A·y + e 这类矩阵乘,本地免费);
- **布尔域**:emp-ag 认证混淆电路(WRK 风格,n 方、不诚实多数、恶意安全),
  负责非线性部分(Decompose、范数比较、拒绝采样判定);
- 两个域之间靠 **edaBits/daBits**(同一值的算术+布尔一致份额对)和
  **A2B 转换**(mask-and-open)过桥。

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

- C++20 编译器、CMake ≥ 3.25、OpenSSL(emp-tool 的依赖)
- **emp 三件套**,均 `cmake --install` 到标准前缀(默认 /usr/local):
  - [emp-tool](https://github.com/emp-toolkit/emp-tool)(基础:NetIO、block、PRG)
  - [emp-ot](https://github.com/emp-toolkit/emp-ot)(OT 扩展)
  - **emp-ag**

每个库的安装方式相同:

```bash
cmake -S . -B build && cmake --build build -j && sudo cmake --install build
```

CMakeLists 里 `find_package(emp-ag REQUIRED)` 会自动带出 emp-tool/emp-ot。

## 跑

```bash
make n=3              # 3 方恶意,编译并本机跑一次 KeyGen + Sign
make n=5              # 5 方(n= 任意值,nP 是模板参数,现编对应二进制)
make n=3 p=44         # 换 ML-DSA-44 参数集(默认 65;还有 p=87)
make n=3 TEST=1       # 测试构建:开出秘密值做范围/电路-明文逐位比对/签名验证
make n=3 TEST=1 TAMPER_C=1   # 对抗测试:party 1 篡改一个 share,MACCheck 应 abort
make clean
```

端口每次随机(避开上一轮 TIME_WAIT);要固定就 `PORT=16400`。

成功输出(绿色为完整签名;红色 `(c, bot, bot)` 是 T=1 下正常的拒绝采样,重跑即可):

```
  pk: rho[0]=9e3779b9, t1=1536 coefficients | tr[0]=...
  signature: (c, z, h) — ||z||_inf=523576 (< 524092), HW(h)=45 (<= 55)
  setup     14.2 ms | comm    0.2 MB
  KeyGen   343.6 ms | comm    7.0 MB | ANDs 17424
  Sign    2021.6 ms | comm  350.3 MB | ANDs 1454736
STATS np=3 param=ML-DSA-65 ...
  KeyGen + Sign complete
```

时间/AND 数/通信量都按阶段拆分(通信量是本方发+收的字节数);`STATS` 行是
同一组数字的机器可读版,供基准脚本 grep。被拒绝的运行 Sign 开销更小
(r₀ 拒绝后跳过 z 电路,AND 数约为接受运行的 60%)。

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
