# ML-DSA Decompose — n 方恶意运行
#
#   make n=3          3 方恶意(emp-ag),编译并跑
#   make n=3 p=44     换 ML-DSA-44 参数集
#   make help         这段说明
#   make clean        删掉 build/
#
# 注:make 不接受 `make --n 3` —— `--` 开头的参数会被 make 当成自己的选项
# (`--n` 还会因为匹配多个内置长选项而报 ambiguous)。用 `n=3`。

# ---- 可调变量 --------------------------------------------------------------
# 注释另起一行:make 会把 `p ?= 65   # 注释` 里 # 前面的空格算进变量值。

# 参与方数。设了就跑 MPC,不设就看 help
n     ?=
# 参数集:44 / 65 / 87
p     ?= 65
# 测试构建:1 会打开秘密值并运行范围/Decompose/adversary 自检
TEST  ?= 0
# 编译 adversary:party 1 在打开前修改一个 c share
TAMPER_C ?= 0
# 1 = slot-restricted (sign_slot.h), 0 = two-round baseline (sign_2round.h)
SLOT  ?= 1
# 环认证的统计安全参数 sigma(论文取 128)
SIGMA ?= 128
# 端口。随机取,避开上一次残留在 TIME_WAIT 里的
PORT  ?= $(shell awk 'BEGIN{srand();print 20000+int(rand()*400)*100}')

BUILD := build
RUNN  := ./runN
# clang-format 不在 PATH(你放的是 .clang-format 配置,不是二进制)。优先用 PATH 里的,
# 否则回退到 Homebrew LLVM。覆盖: make format CLANG_FORMAT=$(xcrun -f clang-format)
CLANG_FORMAT ?= $(shell command -v clang-format 2>/dev/null || echo /opt/homebrew/opt/llvm/bin/clang-format)
PARAM := MLDSA_$(p)
# 默认建 2/3/5,再加上本次 n= 要的那个
PARTIES := 2;3;5$(if $(n),;$(n))

ifdef n
.DEFAULT_GOAL := run
else
.DEFAULT_GOAL := help
endif

.PHONY: help run configure clean format

# ---- 说明 ------------------------------------------------------------------
help:
	@echo 'ML-DSA Decompose'
	@echo ''
	@echo '  make n=3          3 方恶意(emp-ag),编译并跑'
	@echo '  make n=3 p=44     换 ML-DSA-44(默认 65)'
	@echo '  make n=3 TEST=1   测试构建:打开 y/e_w/w/w0 并做自检'
	@echo '  make n=3 TEST=1 TAMPER_C=1   编译 adversary share 篡改测试'
	@echo '  (系数个数 N 硬编码在 decompose.h: constexpr int N)'
	@echo '  make clean        删掉 build/'
	@echo '  make format       clang-format 就地格式化(用 mldsa/.clang-format)'
	@echo ''
	@echo '  端口自动随机(本次会用 $(PORT));要固定就 PORT=16400'

# ---- cmake ------------------------------------------------------------------
# cmake 幂等且快,每次都跑一遍,省得 n= 变了忘记重新 configure。
configure:
	@cmake -S . -B $(BUILD) -DCMAKE_BUILD_TYPE=Release \
	       -DMLDSA_PARTIES="$(PARTIES)" -DMLDSA_PARAM=$(PARAM) \
	       -DMLDSA_TEST=$(if $(filter 1,$(TEST)),ON,OFF) \
	       -DMLDSA_TAMPER_C=$(if $(filter 1,$(TAMPER_C)),ON,OFF) \
	       -DMLDSA_SLOT=$(if $(filter 1,$(SLOT)),ON,OFF) -DMLDSA_SIGMA=$(SIGMA) \
	       -DMLDSA_DRIVER=ON -DMLDSA_TESTS=OFF >/dev/null
	@# p=/TEST= 换配置只改 -D,cmake 的 Makefile 生成器靠 mtime 判断要不要重编,
	@# 而重新生成的 flags.make 可能和上一次的 .o 落在同一秒 —— make 只在依赖
	@# 严格更新时才动手,于是静默沿用旧二进制。配置真变了就删掉 .o。
	@if [ "$$(cat $(BUILD)/.param 2>/dev/null)" != "$(PARAM)-TEST$(TEST)-TAMPER_C$(TAMPER_C)-SLOT$(SLOT)-S$(SIGMA)" ]; then \
	    rm -f $(BUILD)/CMakeFiles/main*.dir/main.cpp.o; \
	    echo $(PARAM)-TEST$(TEST)-TAMPER_C$(TAMPER_C)-SLOT$(SLOT)-S$(SIGMA) > $(BUILD)/.param; \
	fi

# ---- n 方恶意 ---------------------------------------------------------------
run: configure
ifndef n
	$(error 没给参与方数。用法: make n=3)
endif
	@cmake --build $(BUILD) --target main$(n) -j
	@echo '--- $(n) 方恶意 (emp-ag), ML-DSA-$(p), TEST=$(TEST), TAMPER_C=$(TAMPER_C), SLOT=$(SLOT), SIGMA=$(SIGMA), 端口 $(PORT) ---'
	@EMP_PORT=$(PORT) $(RUNN) $(n) ./$(BUILD)/main$(n)

clean:
	rm -rf $(BUILD)

# ---- 格式化 -----------------------------------------------------------------
# 就地格式化本项目所有 C++ 源(用这里的 .clang-format;跳过 build/)。
format:
	@command -v $(CLANG_FORMAT) >/dev/null 2>&1 || test -x "$(CLANG_FORMAT)" || { \
	    echo 'clang-format 未找到: $(CLANG_FORMAT)'; \
	    echo '  装: brew install clang-format   或: make format CLANG_FORMAT=$$(xcrun -f clang-format)'; \
	    exit 1; }
	@find . -type f \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' -o -name '*.cc' \) \
	    -not -path './build/*' -print0 | xargs -0 $(CLANG_FORMAT) -i
	@echo 'clang-format 完成 ($(CLANG_FORMAT))'
