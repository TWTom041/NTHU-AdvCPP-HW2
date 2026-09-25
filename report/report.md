<div class="title" markdown="1">
<div class="t1">Advanced C++ Homework 2</div>
<div class="t2">學號（請填寫）, 姓名（請填寫）</div>
<div class="t2" style="font-size:11pt">Floating Point 行為專案報告 —— 以 Kerr 黑洞測地線模擬引擎為例</div>
</div>

> **報告情境。** 我要為團隊寫一個「黑洞影像」的科學計算程式：在旋轉黑洞（Kerr 時空）附近追蹤數十萬條光線，算出望遠鏡看到的影像，並與天文界實際使用的 [RAPTOR](https://github.com/tbronzwaer/raptor)（Bronzwaer et al. 2018, *A&A* 613, A2；為事件視界望遠鏡 EHT 產生合成影像的 ray tracer 之一）比對。主管擔心「浮點數不安全」。本報告以這個真實問題為主軸：先說明我實作的引擎與驗證方式（第 1 章），再用它以及針對性的小實驗回答 Q0–Q7。**所有程式碼、原始終端輸出與重現腳本都在 repository 中**，每一個數據都可以在 `results/*.txt` 找到出處；截圖是這些輸出的畫面。

**一句話結論：** 浮點數不是「不安全」，而是**有規格、可預測、但不是實數**。只要知道誤差從哪裡來（捨入、抵銷、問題本身的病態性、編譯器/函式庫/硬體的差異），就能量化並控制它 —— 本報告的引擎用 `double` 就能把黑洞陰影邊界算到與解析解相差 10⁻¹²，而 `-ffast-math`、`float16`、跨平台差異等陷阱也都能用實驗重現並說明原因。

[TOC]

<div class="pagebreak"></div>

## Q0. 實驗環境

所有實驗在同一台雲端 Linux 虛擬機上執行（Q5 另外使用 aarch64 交叉編譯 + QEMU 與 musl libc 模擬「不同環境」）。

| 項目 | 規格 |
|---|---|
| CPU | Intel Xeon Processor @ 2.10 GHz（family 6, **model 207 = 0xCF, Emerald Rapids / Raptor Cove 核心**，stepping 2），KVM 虛擬機 4 vCPU（1 thread/core） |
| 指令集 | SSE4.2, AVX2, **FMA**, **AVX-512** (F/DQ/BW/VL/IFMA/VBMI…), **AVX512-FP16**, AVX512-BF16, AMX |
| 快取 | L1d 48 KiB/核, L1i 32 KiB/核, L2 2 MiB/核, L3 260 MiB（共享） |
| 記憶體 | 15 GiB |
| OS / kernel | Ubuntu 24.04.4 LTS, Linux 6.18.44（KVM guest） |
| C 函式庫 | **glibc 2.39**（Ubuntu 2.39-0ubuntu8.7）；Q5 另用 **musl 1.2.4** |
| 編譯器 | **GCC 13.3.0**（Ubuntu 13.3.0-6ubuntu2~24.04.1），`-std=gnu++23`；Q5 用同版本 **aarch64-linux-gnu-gcc 13.3.0** |
| 其他 | QEMU user 8.2.2（執行 aarch64 binary）、Valgrind 3.22.0（Cachegrind）、MPFR 4.2.1（正確捨入參考值）、GSL 2.7.1（RAPTOR 需要）、oneTBB 2021.11（`std::execution::par` 後端）、Python 3.11 / NumPy 2.4 / Matplotlib 3.11、Node 22 + Playwright/Chromium（截圖與 PDF） |
| RAPTOR | github.com/tbronzwaer/raptor，固定 commit `08cb9a2`（2021-05-25），**原始碼未修改** |

**編譯參數**：每個實驗的參數明列在 `Makefile`（因為有些實驗本身就在研究參數），預設為 `g++ -std=gnu++23 -O2`；benchmark 另外標示 `-O3 -march=native`、`-ffast-math` 等。

**量測方法**：`std::chrono::steady_clock`，每項量測重複 3 次取最小值；單執行緒量測時不跑其他工作。雲端 VM 的 turbo 頻率不可控，因此只比較**相對**數字；另外這台 VM **沒有開放硬體效能計數器**（`perf stat` 的 cache 事件回報 `<not supported>`），所以快取命中率用 Valgrind Cachegrind 模擬、再以實測的 working-set sweep 佐證。

**重現**：`make raptor && make all cross && ./scripts/run_all.sh`（約 30 分鐘），`results/` 內即為本報告所有數據；`python3 scripts/make_figures.py && node scripts/screenshots.js && python3 scripts/build_report.py` 產生圖、截圖與本 PDF。

![](shots/q0_env.png){.shot}
<div class="caption">圖 0：實驗環境（lscpu 節錄）</div>

<div class="pagebreak"></div>

## 1. 專案：Kerr 黑洞測地線模擬引擎

### 1.1 物理問題

光在黑洞附近沿著**零測地線 (null geodesic)** 前進。從相機的每個像素往回追蹤一條光線：它可能掉進事件視界（黑色「陰影」）、打到吸積盤（發光）、或逃到無窮遠。旋轉黑洞由 Kerr 度規描述（單位 G = c = M = 1，自旋 a）。

### 1.2 我的實作（`include/kerr/`，header-only C++23）

| 設計 | 說明 |
|---|---|
| **座標** | ingoing Kerr–Schild (t, r, θ, φ)：在事件視界上沒有座標奇異，光線可以平順穿越（Boyer–Lindquist 座標在視界處 g_rr → ∞）。 |
| **方程式** | **Hamiltonian 形式** H = ½ g<sup>ab</sup>(x) p<sub>a</sub>p<sub>b</sub> = 0，ẋ<sup>a</sup> = g<sup>ab</sup>p<sub>b</sub>，ṗ<sub>a</sub> = −∂H/∂x<sup>a</sup>。KS 的逆度規只依賴 (r, θ)，所以 **p<sub>t</sub>（能量）與 p<sub>φ</sub>（角動量）在程式裡是逐位元守恆的**；只剩 H 與 Carter 常數 Q 會漂移，拿來當精度診斷。RAPTOR 則是直接積分二階測地線方程，每步計算 64 個 Christoffel 符號。 |
| **型別** | 所有程式碼以 `template <class T>` 撰寫，可實例化為 `float`、`double`、`long double`、`__float128`、`std::float16_t`，甚至 `std::experimental::simd<double>`（一次推進 8 條光線）。數學函式透過 `kerr::m::sin` 等多載統一（`__float128` 對應 libquadmath 的 `sinq`…）。 |
| **積分器** | RK2（RAPTOR 預設）、RK4、Dormand–Prince 5(4)（自適應步長、FSAL）。步長規則與 RAPTOR 的 `stepsize()` 相同：dλ = h / (\|ṙ\|/r + \|θ̇\|/min(θ, π−θ) + \|φ̇\|)。 |
| **相機** | 與 RAPTOR 的 `initialize_photon()` 相同：像素 → 撞擊參數 (α, β) → (E, L, Q)，在 r = 10⁴ M 由 Boyer–Lindquist 轉到 KS。 |
| **場景** | 幾何薄、光學厚的 Keplerian 吸積盤（r<sub>ISCO</sub> ~ 20 M），觀測強度 I = g⁴ · r⁻³(1 − √(r<sub>in</sub>/r))，g 為重力紅移 + 都卜勒因子。 |

![](fig/render.png)
<div class="caption">圖 1：左：本引擎產生的影像（a = 0.9375、傾角 60°、double、RK4）；中：被黑洞捕獲的光線（黑）與 Bardeen 解析陰影邊界（橘線）；右：每條光線的積分步數 —— 貼近光子環的光線要多繞好幾圈，成本是一般光線的 2–6 倍（這在第 1.5 節的平行化中很重要）。</div>

### 1.3 驗證一：與解析解比對（`src/kerr_validate.cpp`）

遠方觀測者看到的陰影邊界有解析解（Bardeen 1973 的 critical curve，由球面光子軌道給出）。我在 12 個方向上用二分法找「捕獲 ↔ 逃逸」的交界，和以 `__float128` 計算的解析曲線比對：

![](shots/validate.png){.shot}

![](fig/validate.png){style="max-width:80%"}
<div class="caption">圖 2：陰影邊界的最大相對誤差。RK4 的誤差隨步長 h 以 h⁴ 下降（3×10⁻⁹ → 3×10⁻¹¹ → 3×10⁻¹³，每次 h 縮小 3 倍、誤差縮小 ≈ 81 = 3⁴ 倍），證明方程式與積分器都正確；自適應 DP45 的極限則由型別決定：float ≈ 2×10⁻⁷（≈ 2 ε<sub>float</sub>）、double ≈ 10⁻¹²、long double ≈ 10⁻¹⁴、__float128 ≈ 10⁻¹⁵（這裡已經被「二分法輸入參數是 double」限制）。</div>

**這是整份報告最重要的一張圖**：它說明「浮點數誤差」其實是可以**量化**的，而且在這個問題中，演算法的截斷誤差（步長）通常遠大於 double 的捨入誤差。

### 1.4 驗證二：與 RAPTOR 逐條光線比對（`src/raptor_compare.cpp`）

`raptor_bridge/raptor_bridge.c` 把 RAPTOR 的 `.c` 檔（不修改）一起 `#include` 成單一 translation unit，直接呼叫它的 `initialize_photon()`、`stepsize()`、`rk4_step()`／`rk2_step()`（RAPTOR 把全域變數定義在 header 裡，GCC ≥ 10 預設 `-fno-common` 會連結失敗，單一 TU 可以繞過；完整的 RAPTOR 程式則用 `-fcommon` 編譯）。對 100×100 條光線（RAPTOR 範例設定：a = 0.9375、i = 60°、視野 40 M、STEPSIZE = 0.03）分別用 RAPTOR、本引擎（同一種積分器、同一個步長規則）與參考解（本引擎、long double、DP45、rtol 10⁻¹³）追蹤，再把逃逸光線推到同一半徑 r = 1.01×10⁴ 比較方向：

![](shots/raptor_rk4.png){.shot}

| RK4, STEPSIZE 0.03 | RAPTOR | 本引擎 |
|---|---|---|
| 初始條件差異（與 RAPTOR） | — | 2.9×10⁻¹⁴（只差捨入） |
| 捕獲光線數（參考解 486） | 486 | 486 |
| 每步時間（單執行緒） | 574 ns | **211 ns**（2.7× 快：Hamiltonian 形式不需 64 個 Christoffel 符號） |
| 能量 E 漂移 / 角動量 L 漂移 | 9.3×10⁻⁶ / 2.3×10⁻⁵ | **0 / 4.5×10⁻¹⁵**（結構上守恆） |
| 最終方向誤差（中位數 / p99） | 9.4×10⁻⁸ / 6.0×10⁻⁶ rad | **8.5×10⁻⁹ / 2.5×10⁻⁸ rad** |

* STEPSIZE 0.01 時兩者都以 h⁴ 收斂（誤差 ≈ 1/81），確認兩個程式解的是同一個問題。
* **RAPTOR 的預設積分器其實是 RK2**（`parameters.h`：`int_method (RK2)`）。在同一步長下 RAPTOR RK2 捕獲 499 條光線（多了 13 條錯誤分類），方向誤差中位數 4×10⁻³ rad、最大 2.5 rad；本引擎的 RK2 為 2.6×10⁻⁴ / 1.2×10⁻²。**換句話說，對這個問題而言誤差的主要來源是積分方法與步長，不是 double 的精度。**
* 最大誤差都出現在陰影邊緣附近 —— 那些光線繞光子環好幾圈，是混沌的（見 Q6）。

![](shots/raptor_rk2.png){.shot}

### 1.5 平行化技術比較（`src/parallel_bench.cpp`）

同一張 256×256 影像（double、RK4、含吸積盤），以不同平行技術渲染，並檢查結果與序列版逐像素一致（mismatch 欄）。機器有 4 核。

![](shots/parallel_double.png){.shot}

![](fig/parallel.png){style="max-width:85%"}
<div class="caption">圖 3：相對於單執行緒 scalar 版本的加速（藍：double，橘：float）。虛線為 4 核的理想值。</div>

觀察與解釋：

1. **負載不平衡決定了排程策略的好壞。** 每條光線的步數從 219 到 1366 不等（圖 1 右），而且昂貴的光線集中在影像中央。把影像切成 4 個連續區塊（`std::thread static block`、`OpenMP schedule(static)`）時，拿到中央區塊的執行緒最慢，imbalance（最忙執行緒時間 / 平均）= 1.36，效率只有 72–73%。`schedule(guided)` 一開始就分大塊，同樣吃虧（72%）。
2. **交錯分配（cyclic rows、`schedule(static,1)`）或動態排程（atomic 工作計數器、`schedule(dynamic,16)`、TBB work stealing 的 `std::for_each(std::execution::par, …)`）都能達到 94–100% 效率（量測雜訊使少數數字略超過 100%）**：前者把昂貴區域平均分給每個執行緒，後者讓先做完的執行緒去拿新工作。
3. **SIMD 是另一個維度**：`std::experimental::simd<double>`（AVX-512，8 lanes）把 8 條相鄰光線綁在一起、以遮罩 (mask) 凍結已結束的 lane。單執行緒就有 **4.6×**（lane 使用率 97%：相鄰像素的步數很接近），再加上 OpenMP 4 核達到 **18.6×**，而且結果與 scalar 版本的差異只有 5×10⁻¹⁴（向量版 `sin`/`cos` 與 glibc 的捨入不同）。
4. **float vs double（Q1 的伏筆）**：scalar 程式中 float 並沒有比較快（這個 `-march=native` build 中 float 甚至慢 27%，見 Q1.4）；但在 SIMD 版本中 float 一個暫存器可以放 16 個 lane：float SIMD 單執行緒是 scalar float 的 **15×**、是 double SIMD 的 2.6×，加上 4 核達到 **60×**，捕獲/逃逸分類完全相同，強度差異 1.7×10⁻⁴（float 的精度）。
5. **False sharing**：4 個執行緒各自遞增「自己的」計數器，但計數器相鄰擺放、落在同一條 64-byte cache line 上：每次寫入都要搶這條 line 的所有權（MESI 協定），**慢了 10 倍**（72 ns vs 7.2 ns 一次遞增）；以 `alignas(64)` 讓每個計數器獨佔一條 line 就解決了。

<div class="pagebreak"></div>

## Q1. Benchmark：各種浮點型別的速度與精度

### Q1.1 型別與格式（IEEE 754）

一個二進位浮點數 = (−1)<sup>s</sup> × 1.f × 2<sup>e−bias</sup>。「精度」由尾數位數 p 決定（machine epsilon ε = 2<sup>1−p</sup>），「範圍」由指數位數決定。

| C++ 型別（x86-64 Linux） | 格式 | 位元（符號/指數/尾數） | p | ε | 最大值 | 十進位有效位數 | 硬體 |
|---|---|---|---|---|---|---|---|
| `std::float16_t` / `_Float16` | binary16 | 1/5/10 | 11 | 9.8×10⁻⁴ | 65504 | ~3 | AVX512-FP16（本機有）；否則只有轉換指令 F16C |
| `float` | binary32 | 1/8/23 | 24 | 1.2×10⁻⁷ | 3.4×10³⁸ | ~7 | SSE/AVX |
| `double` | binary64 | 1/11/52 | 53 | 2.2×10⁻¹⁶ | 1.8×10³⁰⁸ | ~16 | SSE2/AVX |
| `long double` | x87 80-bit extended | 1/15/64（顯式整數位） | 64 | 1.1×10⁻¹⁹ | 1.2×10⁴⁹³² | ~19 | x87 FPU（佔 16 bytes） |
| `__float128` | binary128 | 1/15/112 | 113 | 1.9×10⁻³⁴ | 1.2×10⁴⁹³² | ~34 | **無硬體**，libgcc 軟體模擬 |

（注意：`long double` 的格式由 ABI 決定 —— x86-64 Linux 是 80-bit、aarch64 Linux 是 binary128、MSVC 是 64-bit，這在 Q5 會成為「結果不同」的原因。）

### Q1.2 浮點數的硬體實作

**(1) 資料路徑。** 現代 x86-64 的 `float`/`double` 運算在 SSE/AVX 暫存器（xmm/ymm/zmm）上進行：scalar 指令（`addss`/`addsd`）只用最低的 lane，packed 指令（`vaddps`/`vaddpd`）一次處理整個暫存器 —— 512-bit 暫存器可放 **16 個 float 或 8 個 double**。本機的 Raptor Cove 核心有兩個 FMA 單元（port 0/5）及兩個「快速加法器」（port 1/5）：加法延遲約 2 cycles、乘法/FMA 約 4 cycles，**float 與 double 的延遲相同**，因為同一組電路以相同級數的 pipeline 處理兩種寬度。

**(2) 捨入。** 加法器先對齊指數、相加得到比 p 位更多的結果，再依 guard / round / sticky 三個額外位元以及 MXCSR 暫存器中的 rounding mode（預設 round-to-nearest-even）捨入回 p 位，並在 MXCSR 設定 sticky 的例外旗標（inexact、underflow、overflow、invalid、divide-by-zero）。**IEEE 754 要求 + − × ÷ √ 與 FMA 都是「正確捨入」**：結果等於把精確值捨入一次 —— 所以基本運算是完全可預測、可重現的。FMA（a×b+c 只捨入一次）是最重要的新增運算。

**(3) 除法與開根號** 使用獨立的迭代式單元（radix-16 SRT/digit recurrence），延遲比乘法長很多且只能部分 pipeline；float 需要的位數少，所以**除法/開根號是 float 比 double 快的少數 scalar 運算**（見下表 div：4.2 vs 5.2 ns）。

**(4) Subnormal（非正規數）** 在很多微架構上需要 microcode assist，一次運算可能慢數十倍（Q6 實測慢 17 倍）；`-ffast-math` 的 FTZ/DAZ 就是為了避開它（Q3）。

**(5) x87 與 `long double`。** x87 是 1980 年的 8087 協處理器架構：8 個 80-bit 暫存器組成**堆疊**，沒有 SIMD，每次載入/存回 80-bit 值都需要 `fld`/`fstp` tbyte（較慢），超越函數（`fsin`）是 microcode。GCC 對 `long double` 產生 `fadd`/`faddp` 等 x87 指令（見截圖）。

**(6) `__float128`：沒有硬體。** x86 沒有 binary128 運算單元，GCC 把每個 `+`、`×` 編譯成對 libgcc 的函式呼叫（`__addtf3`、`__multf3`…），在函式裡用 64-bit 整數運算拆解指數與尾數、處理捨入 —— 這就是 soft-float。

**(7) `float16`。** 在基本指令集（x86-64-v1）下沒有 half 運算，GCC 為每個運算呼叫 `__extendhfsf2`（轉 float）→ float 運算 → `__truncsfhf2`（轉回 half）；本機支援 **AVX512-FP16**，用 `-march=native` 時 GCC 直接產生 `vaddsh`/`vaddph`。

![](shots/q1_asm.png){.shot}
<div class="caption">圖 4：反組譯證據（上到下）：float16 在 -O2 呼叫轉換函式、在 -march=native 使用 AVX512-FP16 指令；long double 使用 x87 `fadd`；__float128 呼叫 libgcc 的 `__addtf3`/`__multf3`。</div>

### Q1.3 Benchmark 1：基本運算（`experiments/q1_microbench.cpp`）

* **latency**：一條相依鏈 `x = x op c`，下一個運算必須等上一個結束 —— 量到的是**單一運算的延遲**（ODE 積分這類有資料相依的程式看到的就是它）。
* **throughput**：`c[i] = a[i] op b[i]`，2048 元素陣列（在 L1 內），運算彼此獨立，編譯器可自動向量化 —— 量到的是**最佳情況下每個元素的成本**。
* 各編譯兩次：`-O2`（x86-64 基本指令集，SSE2，128-bit 向量）與 `-O3 -march=native`（AVX-512）。

![](shots/q1_micro.png){.shot}

![](fig/q1_micro.png)
<div class="caption">圖 5：ns/operation（對數座標，越低越好）。</div>

**解釋：**

1. **scalar 的 float 與 double 延遲幾乎一樣**（add 0.71 vs 0.72 ns、mul 1.6 vs 1.5 ns）：同一組硬體（上一節 (1)）。只有 div/sqrt 因迭代位數不同，float 稍快。
2. **向量化時 float 是 double 的 2 倍**：`-O2` 的 add 0.136 vs 0.202 ns、`-march=native` 0.057 vs 0.116 ns —— 一個暫存器裝的 float 數量是 double 的兩倍；sqrt 在兩種 build 都是 1.07 vs 2.16 ns。**float 的速度優勢來自 SIMD 寬度與記憶體頻寬（Q1.5），不是單一運算比較快。**
3. **`long double`**：加乘延遲接近 double（x87 本身不慢），但無法向量化（throughput 3.6 ns，比 double 慢 30 倍），而且 `sinl` 的實作較貴（35 ns，見 Q2）。
4. **`__float128`**：每個運算都是函式呼叫 + 整數模擬，加法 18 ns、除法 28 ns、**sqrt 135 ns、sin 379 ns**（比 double 慢 25–60 倍）。
5. **`float16_t`**：在 -O2 每個運算 13 ns（兩次轉換函式呼叫！），比 double 慢 20 倍；開 `-march=native` 後有硬體支援，向量加法 0.035 ns/元素（**一個 zmm 放 32 個 half**），是所有型別中最快的 —— 但精度只有 3 位數（下一節）。

### Q1.4 Benchmark 2：真實工作量 —— 同一個 Kerr 光線追跡（`src/kerr_precision.cpp`）

48×48 條光線（含吸積盤，RK4），單執行緒。**參考解是同一個演算法用 `__float128` 執行**，所以誤差欄位是**純捨入誤差**（截斷誤差相同而抵銷）。

![](shots/q1_kerr.png){.shot}

| 型別 | ns/step | 相對 double | 影像強度誤差（中位數 / 最大） | 光線方向誤差（中位數） | 評語 |
|---|---|---|---|---|---|
| float | 200 | 0.90× 時間 | 9.0×10⁻⁷ / 7.1×10⁻⁵ | 1.3×10⁻⁶ rad | scalar 程式幾乎不快；誤差 ~ ε<sub>float</sub> × 放大 |
| double | 222 | 1 | 1.8×10⁻¹⁵ / 9.2×10⁻¹⁴ | 2.6×10⁻¹⁵ rad | 捨入誤差遠小於截斷誤差（~10⁻⁸） |
| long double | 1508 | 6.8× 慢 | 8.6×10⁻¹⁹ | 1.3×10⁻¹⁸ rad | x87，無 SIMD；`sinl`/`cosl` 慢 |
| __float128 | 19705 | 89× 慢 | 0（參考解） | 0 | 軟體模擬 |
| float16_t | — | — | 全部 2304 條光線失敗 | — | **溢位**：相機設定中 (r²+a²)² ≈ 10¹⁶ ≫ 65504 → ∞ → NaN |

解釋：

* **為什麼 float 在這裡沒有變快？** 光線追跡是 scalar、相依鏈很長（每一步依賴上一步）、每次 RHS 呼叫 `sin`/`cos` 的程式；如 Q1.3 所示，scalar float 與 double 的延遲相同，`sinf` 只比 `sin` 快一點。記憶體用量極小（見 Q1.5 的 Cachegrind），也沒有頻寬優勢。唯有像第 1.5 節那樣**把資料排成 SIMD**，float 才有 2 倍的空間。
* **編譯器也會讓 float 變慢**：用 `-march=native` 編譯時，GCC 的 SLP 自動向量化把 8 個元素的 RK4 狀態陣列打包成向量，產生的 shuffle 反而拖慢程式：float 從 223 變成 276 ns/step（double 225 → 225）；加上 `-fno-tree-vectorize` 後 float 220、double 179 ns/step。「float 比 double 快」並不是硬體的保證，而是取決於資料排列與編譯器。
* **H 漂移**（max\|2H\|）：double、long double、`__float128` 都是 3.4×10⁻⁵ —— 這是 RK4 的截斷誤差，和型別無關；float 為 1.7×10⁻³，捨入誤差已經超過截斷誤差。**這告訴我們 double 對這個問題是「精度剛好過剩」：加更多位數（long double、quad）完全沒有幫助，卻付出 7–90 倍的時間。**
* **float16 的失敗是「範圍」問題不是「精度」問題**：中間值 (r²+a²)² 在 r = 10⁴（甚至 r = 100）時就超過 65504。選擇型別時除了有效位數，還要檢查中間運算的動態範圍。

### Q1.5 記憶體與快取：float 真正的優勢在哪裡

**(a) Working-set sweep（實機量測，`experiments/cache_sweep.cpp`）**：對大小 4 KiB–1 GiB 的陣列做加總（float vs double），以及隨機指標追逐（測 load-to-use 延遲）。

![](fig/cache.png)
<div class="caption">圖 6：左：每 ns 加總的元素數；右：相依載入延遲。虛線為 L1d/L2/L3 容量。</div>

* 指標追逐清楚顯示階層：**L1 ≈ 2 ns、L2 ≈ 6–9 ns、L3 ≈ 33–46 ns（2–8 MiB）**；超過 ~16 MiB 後延遲逐步升到 **150–215 ns**（VM 實際可用的 L3 份額與 4 KiB 分頁的 TLB miss 都會貢獻）。
* 資料在 L1/L2 內時，加總受限於運算，float（16 lanes）每 ns 處理約 **1.5–1.9 倍**於 double 的元素；超出快取後受限於記憶體頻寬（兩者都約 13 GB/s），float 因為**只搬一半的位元組**而處理 **2 倍**的元素（3.3 vs 1.7 Gelem/s）。

**(b) Cache 命中率（Valgrind Cachegrind 模擬，`experiments/cache_patterns.cpp`）**：沿列 (row-major, stride 1) 與沿行 (column-major, stride N) 加總 2048×2048 矩陣，以及光線追跡程式本身。

![](shots/cachegrind.png){.shot}

* 沿列走訪時，每條 64-byte cache line 含 16 個 float／8 個 double，**D1 讀取 miss rate 5.8% 與 11.4%（≈ 1/16 與 1/8）**（float 命中率較高，這就是「float 省頻寬」的微觀原因）；沿行走訪每次都碰新的一條 line，D1 讀取 miss rate **89%**，實測（4096×4096）慢 9–10 倍。
* **光線追跡程式的 D1 miss rate 接近 0**：每條光線的狀態只有 8 個數、幾十個暫存值，全部在暫存器與 L1 中 —— 這是一個 **compute-bound** 程式，這也解釋了為何 Q1.4 中 float 沒有頻寬優勢。

### Q1.6 小結：如何選型別

1. **預設用 `double`**：scalar 速度與 float 相同，精度多 9 位數，範圍 10³⁰⁸ 幾乎不會溢位。
2. **`float` 只在「可向量化或受記憶體頻寬限制」且誤差預算 ≥ 10⁻⁶ 時才划算**（影像、圖學、ML、大陣列）。
3. **`long double`/`__float128` 用來驗證**（例如本報告的參考解），不要當作主要型別：慢 7–90 倍，而且 `long double` 在不同平台是不同格式。
4. **`float16`/`bfloat16` 只適合儲存或容錯的 ML 計算**：3 位有效數字、最大值 65504。
5. 比「換更長的型別」更有效的通常是**換更好的演算法**（RK4 vs RK2、守恆形式的方程式、避免抵銷的公式）。

<div class="pagebreak"></div>

## Q2. Library Implementation：`sin(x)`

### Q2.1 規格

* **C++**：`<cmath>` 的 `std::sin` 規格引用 C 標準 7.12.4.6：「sin 函數計算 x（以弧度為單位）的正弦」。C++ 提供 `float sin(float)`、`double sin(double)`、`long double sin(long double)`（以及 `sinf`、`sinl`），整數引數當作 `double`；C++23 增加 extended floating-point 型別（`std::float16_t` 等）的多載，C++26 起為 `constexpr`（P1383）。
* **精度**：C 標準 5.2.4.2.2 說 `<math.h>` 函數的精度是 **implementation-defined** —— 標準**不要求**正確捨入，也不要求不同實作結果相同。IEEE 754-2019 §9.2 只「建議 (should)」sin 正確捨入。
* **特殊值**（C Annex F.10.1.6，glibc 遵守）：`sin(±0) = ±0`（保留符號）；`sin(±∞)` 回傳 NaN 並觸發 **FE_INVALID**，若 `math_errhandling & MATH_ERRNO` 則 `errno = EDOM`；`sin(NaN) = NaN`；很小的非正規數回傳 x 本身並觸發 underflow。

![](shots/q2_special.png){.shot}

### Q2.2 glibc 2.39 的實作（`sysdeps/ieee754/dbl-64/s_sin.c`，源自 IBM Accurate Mathematical Library）

x86-64 的 glibc **不使用 x87 的 `fsin` 指令**，完全以 SSE2 的 double 運算實作。依 |x| 的高 32 位元分成 5 個區間：

| 區間 | 做法 |
|---|---|
| \|x\| < 2⁻²⁶ | 直接回傳 x（sin x = x − x³/6，x³/6 小於半個 ulp） |
| 2⁻²⁶ ≤ \|x\| < 0.855469 | `do_sin(x, 0)`：若 \|x\| < 0.126 用 Taylor 多項式；否則**查表**：把 x 拆成 x = x<sub>i</sub> + dx，x<sub>i</sub> 是 1/128 的倍數（技巧：`u.x = big + |x|`，big = 1.5·2⁴⁵，加法的捨入自動把 x 四捨五入到 2⁻⁷ 的倍數，低位元就是表格索引）。表 `__sincostab` 有 440 個 double，存 sin(x<sub>i</sub>)、cos(x<sub>i</sub>) 以及它們的「尾數修正項」(double-double)。再用 sin(x<sub>i</sub>+dx) = sin x<sub>i</sub> cos dx + cos x<sub>i</sub> sin dx，dx 很小，sin dx、cos dx 只要 3 項多項式。 |
| 0.855469 ≤ \|x\| < 2.426265 | 用 sin x = cos(π/2 − x)：π/2 以兩個 double（hp0 + hp1，共 106 位元）表示，t = hp0 − \|x\| 精確，再呼叫 `do_cos(t, hp1)` |
| 2.426265 ≤ \|x\| < 105414350 | **Cody–Waite 範圍縮減** `reduce_sincos`：n = round(x · 2/π)（技巧：加上 toint = 1.5·2⁵² 讓 FPU 幫忙取整）；y = x − n·(π/2)，其中 π/2 拆成 4 段（mp1、mp2、pp3、pp4），前幾段只有 20–30 個有效位元，使 n·mp1 在 n < 2²⁷ 時**精確**（無捨入），總共約 136 位元的 π/2，得到 double-double 的 (a, da)；再依 n mod 4 選 ±sin/±cos |
| \|x\| ≥ 105414350 | **`__branred`（Payne–Hanek 式）**：x 可能大到 10³⁰⁸，x·(2/π) 需要 2/π 的第 1000 多位才能得到正確的小數部分。表 `toverp` 存 2/π 的二進位展開（每段 24 位元），只取和 x 的指數對齊的那幾段相乘累加，得到 x mod π/2 |
| Inf / NaN | 回傳 x/x（產生 NaN 與 FE_INVALID），Inf 時設 errno = EDOM |

其他實作細節：

* **捨入模式**：函式開頭 `SET_RESTORE_ROUND_53BIT(FE_TONEAREST)` —— 若呼叫者用 `fesetround` 改了捨入模式，函式內暫時切回 round-to-nearest，演算法的誤差分析才成立。
* **錯誤界**：原始碼註解標示各區間 max ULP 0.501–0.548（「~0.55 ULP」）—— **不是正確捨入**（正確捨入 ≤ 0.5 ULP）。
* **硬體 dispatch（IFUNC）**：`sysdeps/x86_64/fpu/multiarch/s_sin.c` 把同一份原始碼編兩次：`__sin_sse2` 與以 `-mfma -mavx2` 編譯的 `__sin_fma`，執行期由 IFUNC resolver 依 CPU 選一個。以 `experiments/q2_ifunc.c` 印出 `sin` 解析後的位址並反組譯：本機預設選到含 36 個 `vfmadd*` 指令的版本；用 `GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA` 遮掉 FMA 後選到另一個位址、0 個 FMA 指令的版本（Q5 會看到這讓結果不同）。

![](shots/q2_ifunc.png){.shot}

### Q2.3 實驗：精度與速度（`experiments/q2_sin.cpp`，以 MPFR 320 位元為正確值）

我另外寫了 60 行的 `mini_sin`（fdlibm 1993 的演算法：3 段 33 位元的 π/2 做 Cody–Waite + 6 次 minimax 多項式），以及直接呼叫 x87 `fsin` 指令，與 glibc 比較。每個區間 20 萬個隨機輸入：

![](shots/q2_sin.png){.shot}

1. **glibc `sin` 在所有區間最大誤差 ≤ 0.515 ULP，99.9% 正確捨入**，包括 x 到 10³⁰⁰ —— 這要歸功於 136 位元的 π/2 與 Payne–Hanek。
2. **x87 `fsin` 在 x = π 誤差 1.6×10¹¹ ULP**：Intel 的 `fsin` 內部只用 66 位元的 π，x 接近 π 的倍數時 x − kπ 的抵銷把 π 的誤差放大 —— 這是 glibc 從 2.28 起完全不使用它的原因（Intel 在 2014 年修正了手冊中的精度宣稱）。
3. **`mini_sin` 在 \|x\| < π/4 附近誤差約 1 ULP**（沒有查表與 double-double 修正），且 **\|x\| 超過約 1.6×10⁶ 後完全錯誤**（1.1×10¹² ULP；x = 10²² 時誤差 3.6×10⁸⁷ ULP）：33 位元一段的 π/2 只保證 n < 2²⁰ 時 n·(π/2 片段) 精確。這說明**範圍縮減才是 sin 最困難的部分**，也說明為什麼 glibc 需要 5 個區間。
4. **速度**：glibc `sin` 在 \|x\| < 2.43 約 6.4 ns，Cody–Waite 區間 23 ns，Payne–Hanek 區間 71 ns（前兩者用了 IFUNC 選出的 FMA 版本）；x87 `fsin` 26–48 ns 且不準；`sinl`（long double）17–179 ns。
5. `sinf` 的最大誤差 0.56 ULP、`sinq`（libquadmath）1.21 ULP、`sinl` 1.40 ULP —— 同一個函數在不同型別的實作精度也不同。

<div class="pagebreak"></div>

## Q3. IEEE 754 與 `-ffast-math`

### Q3.1 它做了什麼

`-ffast-math` 告訴 GCC：「**我不在乎 IEEE 754 的語意，請把浮點數當實數最佳化**」。它等同於開啟下列選項，並定義巨集 `__FAST_MATH__`：

| 選項 | 允許編譯器做的事 | 破壞了什麼 |
|---|---|---|
| `-fno-math-errno` | `sqrt` 等不設定 errno → 可以直接用 `sqrtsd` 指令、可向量化 libm 呼叫 | errno 錯誤回報 |
| `-funsafe-math-optimizations` = `-fassociative-math` + `-freciprocal-math` + `-fno-signed-zeros` + `-fno-trapping-math` | (a+b)+c → a+(b+c)（可向量化歸約）；x/y → x·(1/y)；x+0.0 → x；假設不會有 trap | 結合律/捨入順序、正確捨入的除法、−0 的符號 |
| `-ffinite-math-only` | 假設沒有 NaN 與 ±Inf | `isnan()`/`isinf()` 可被直接判定為 false |
| `-fno-rounding-math`, `-fno-signaling-nans` | 假設預設捨入模式、無 signaling NaN | （本來就是 GCC 預設） |
| `-fcx-limited-range` | 複數除法不做縮放與 NaN 修正 | 複數除法的溢位與 Inf 處理 |
| `-fexcess-precision=fast` | 允許中間結果保留在較寬的暫存器 | 在 x87/float16 上結果可能依暫存器分配而變 |
| 連結 `crtfastmath.o` | **程式啟動時設定 MXCSR 的 FTZ（flush-to-zero）與 DAZ（denormals-are-zero）** | 非正規數全部變成 0 —— **影響整個行程**，包括沒用 `-ffast-math` 編譯的函式庫 |

（`-Ofast` = `-O3 -ffast-math` + 其他；**RAPTOR 的 makefile 就是用 `-Ofast`**。）

### Q3.2 範例：只因為 `-ffast-math` 而改變的結果（`experiments/q3_fastmath.cpp`）

同一份原始碼以 `g++ -O3` 與 `g++ -O3 -ffast-math` 各編譯一次；所有輸入來自 `volatile`，避免常數摺疊：

![](shots/q3_strict.png){.shot}
![](shots/q3_fast.png){.shot}

| # | 實驗 | `-O3` | `-O3 -ffast-math` | 原因 |
|---|---|---|---|---|
| 1 | 陣列中有 0.0/0.0 產生的 NaN，`std::isnan` 檢查 | true | **false** | `-ffinite-math-only`：isnan 被最佳化成常數 false —— **錯誤檢查整段消失** |
| 1b | `std::isinf(1e308 * 10)` | true | **false** | 同上 |
| 2 | 10⁷ 個 0.1f 加總：naive / **Kahan** | 1087937 / **1000000** | 975599 / **975599** | `-fassociative-math` 認為 Kahan 的補償項 c = (t − s) − y 代數上恆為 0，把它刪掉，**Kahan 退化成普通加總**；naive 加總也因為向量化改變順序而得到不同的錯誤結果 |
| 3 | 2²⁰ 項不同數量級的內積 | …8898535 | …8382154 | 歸約被重新結合成 SIMD 的部分和，捨入順序不同 |
| 4 | `1/(-0.0 + 0.0)` | +inf | **−inf** | `-fno-signed-zeros`：x + 0.0 → x，但 −0 + +0 在 IEEE 是 +0 |
| 5 | 10⁶ 個 x/3 | 0 個與正確捨入不同 | **333328 個不同** | `-freciprocal-math`：x/3 → x·(1/3)，兩次捨入，約 1/3 的結果差 1 ulp |
| 6 | 1e-308 × 1e-5（應為非正規數 1e-313） | 1e-313 | **0** | `crtfastmath.o` 設了 FTZ/DAZ（MXCSR 0x1f80 → 0x9fc0） |
| 7 | `sqrt(-1)` 的 errno | EDOM | **0** | `-fno-math-errno` |
| 8 | 內積速度 | 0.83 ns/元素 | 0.59 ns/元素 | 允許向量化歸約（1.4× 快） |

**最典型的單一例子**是 #2：Kahan summation 是教科書上「修正浮點誤差」的演算法，而 `-ffast-math` 會**默默把它最佳化掉** —— 因為在實數的世界裡它什麼都沒做。#1 則是最危險的：程式中防呆的 NaN 檢查會消失，錯誤在下游才以奇怪的形式出現。

### Q3.3 真實程式：Kerr 引擎與 RAPTOR

![](shots/q3_kerr.png){.shot}

* **本引擎**（192×192，含吸積盤）：double 版本有 **40% 的像素逐位元改變**（最大相對差 2.4×10⁻¹³，總通量在第 16 位不同）；float 版本 43% 的像素改變（最大相對 1.3×10⁻⁴、總通量差 10⁻⁷）；捕獲/逃逸分類都沒有改變。速度差異 < 2% —— 這個程式的瓶頸是相依鏈與 `sin`/`cos`，`-ffast-math` 幾乎沒有好處，卻讓結果不可重現。
* **RAPTOR 預設的 `-Ofast`** vs 符合 IEEE 的 `-O3`（完整 GRMHD 影像 500×500，dump040，100 GHz）：總通量印出的 7 位數相同，但有 30 個像素在輸出的 6 位有效數字上不同，**最大相對差異 3.7%**（都在少數貼近光子環或低強度的像素）；執行時間 69 s vs 72 s，只快 4%。以光線層級比較（`raptor_compare_ofast`），`-Ofast` 下 RAPTOR 的測地線統計數字與 `-O2` 完全相同。

![](fig/raptor.png){style="max-width:80%"}
<div class="caption">圖 7：左：RAPTOR 產生的 GRMHD 影像（-O3）；右：-Ofast 改變的像素位置與相對大小。</div>

![](fig/fastmath.png){style="max-width:50%"}
<div class="caption">圖 8：本引擎 float 版本在 -ffast-math 下改變的像素（相對差異，對數色階）。</div>

### Q3.4 建議

* **不要對整個專案開 `-ffast-math`**；它會關掉 NaN 檢查、破壞補償演算法，並透過 FTZ/DAZ 影響同一行程中的其他程式碼。
* 需要速度時**只開真正需要的子選項**：`-fno-math-errno`（幾乎總是安全，讓 `sqrt` 內聯）、`-fno-trapping-math`；需要向量化歸約時自己寫多個部分和（如 `experiments/cache_sweep.cpp`），而不是讓編譯器重新結合。
* 如果一定要用，對關鍵函式以 `__attribute__((optimize("no-fast-math")))` 或把它放到獨立的編譯單元。

<div class="pagebreak"></div>

## Q4. 浮點數輸出：`printf("%.9f")` 是四捨五入嗎？

**結論：這個敘述不精確。** `printf("%.9f", x)` 與 `std::cout << std::fixed << std::setprecision(9) << x`（libstdc++ 的 `num_put` 依 C++ 標準 [facet.num.put.virtuals] 以 printf 轉換規則實作）做的是：**把 x 所代表的「精確二進位值」依「目前的捨入模式」捨入到小數第 9 位**。C 標準 7.23.6.1 只說 "the value is rounded"，Annex F（IEC 60559）要求二進位→十進位轉換依目前捨入方向正確捨入；glibc 與 musl 都是精確實作。和「把我寫的十進位數字四捨五入」有四個差異：

1. **被捨入的是二進位值，不是你寫的十進位字面值。** `0.1234567895` 存成 0.12345678949999999707…，所以印出 `0.123456789`（四捨五入預期 …790）；`1.0000000015` 存成 1.00000000149999990…，印出 `1.000000001`（預期 …002）。反過來 `1.0000000005` 剛好存成比 …5 大一點的值，所以結果「看起來」對 —— 對不對取決於運氣。
2. **剛好在中間（tie）時是「四捨六入五成雙」(round-half-to-even)，不是四捨五入。** 1/1024 = 0.0009765625 **可以被 double 精確表示**，是真正的 tie：`%.9f` 印出 `0.000976562`（2 是偶數，捨去），四捨五入應為 `0.000976563`；5/1024 = 0.0048828125 → `0.004882812`。同理 `%.0f` 對 0.5、1.5、2.5、3.5 印出 0、2、2、4。
3. **結果依動態捨入模式 (`fesetround`) 改變。** 同一個 0.1，在 `FE_UPWARD` 下印出 `0.100000001`；0.3 在 `FE_DOWNWARD`/`FE_TOWARDZERO` 下印出 `0.299999999`；−0.3 在 `FE_UPWARD` 下印出 `-0.299999999`。iostream 也一樣。
4. **負數捨入到 0 時保留負號**：−1e-10 印出 `-0.000000000`，而不是 `0.000000000`。

另外 `float` 會先被提升為 double 再印出它的精確值：`0.1f` 印出 `0.100000001`（0.1f = 0.100000001490116…），而 `16777217.0f` 根本無法表示，印出 `16777216.000000000`。

![](shots/q4_glibc.png){.shot}
<div class="caption">圖 9：glibc 的輸出（每列最後為以 %.40f 印出的精確二進位值；iostream 欄「same」表示與 printf 相同）。musl 1.2.4 的結果完全相同（results/q4_printf.txt）。</div>

**若真的需要「十進位四捨五入」**（例如金額），應使用十進位型別或整數（以「分」為單位），或先用 `std::to_chars` 取得最短可回讀的十進位字串再以十進位規則捨入 —— 而不是依賴 `printf`。

<div class="pagebreak"></div>

## Q5.（加分）相同程式碼、相同編譯器版本、相同參數、無 UB，結果卻不同

`experiments/q5_env.c` 用**同一份原始碼、同為 GCC 13.3.0、同樣只加 `-O2`** 編出三個版本：x86-64 + glibc 2.39、x86-64 + musl 1.2.4（`musl-gcc` 只是包裝同一個 gcc）、aarch64 + glibc 2.39（`aarch64-linux-gnu-gcc`，以 QEMU 執行）。程式沒有未定義行為，所有 libm 的輸入以 `volatile` 強制捨入，保證三個平台餵給函式庫的位元完全相同。

![](shots/q5_platforms.png){.shot}

### 例一：硬體（ISA）不同 → FMA 自動融合（contraction）

`double mul_sub(a,b,c) { return a*b - c; }`，`0.1*10 - 1`：x86-64 印出 **0**、aarch64 印出 **5.55×10⁻¹⁷**。原因：GCC 在 GNU 模式預設 `-ffp-contract=fast`，**只要目標指令集有 FMA 就把 a*b−c 融合成一條 FMA**（只捨入一次）。aarch64 的基本指令集就有 FMA（反組譯為 `fnmsub d0,d0,d1,d2`），x86-64 的基本指令集（SSE2）沒有，所以是 `mulsd` + `subsd`（兩次捨入：0.1×10 先捨入成剛好 1.0）。C/C++ 標準允許這種 contraction（C 6.5p8、`FP_CONTRACT`），所以不是 UB。後果不只最後一位：二次方程式重根的判別式 b²−4ac 在 x86-64 為 0（一個重根），在 aarch64 為 8.9×10⁻¹⁸，開根號得到 3×10⁻⁹ —— **重根變成兩個不同的根**。

同樣的原因，**我的 Kerr 引擎**在 aarch64 上（同樣 `-O2`）有 4054/9216 = 44% 的像素與 x86-64 逐位元不同（最大相對差 5.5×10⁻¹⁴，分類結果相同）。

### 例二：ABI 不同 → `long double` 是不同的型別

`(1.0L + 1e-17L) - 1.0L`：x86-64 為 9.97×10⁻¹⁸（80-bit，64 位尾數）、aarch64 為 1×10⁻¹⁷（binary128，113 位尾數）；十個 0.1L 相加與 1.0L 的差：1.08×10⁻¹⁹ vs −9.6×10⁻³⁵。（MSVC 上 `long double` 就是 double，又會是第三種答案。）

### 例三：作業系統 / C 函式庫不同 → libm 不同

C 標準不要求 `sin`、`exp` 正確捨入（Q2），各家實作的演算法不同。同樣是 x86-64、相同輸入，glibc 與 musl 的結果：

![](shots/q5_diffs.png){.shot}

musl 的 `sin`/`cos`/`tan` 源自 fdlibm，約 3% 的結果與 glibc 差 1 ulp；`atan2` 18%、`cbrt` 50%（最大 3 ulp）。有趣的是 musl 的 `exp`/`log`/`pow` 與 glibc 的差異，和 glibc「FMA 版 vs 非 FMA 版」的差異**完全相同**（611/1/726 個）—— 兩者用的是同一份演算法（ARM optimized-routines），差別只在有沒有 FMA。aarch64 的 glibc 與 x86-64 的 glibc（FMA 版）在 `sin`/`exp`/`pow` 等完全一致，但 `cbrt`（19%，最大 4 ulp）與 `lgamma`（23%）不同 —— 這兩個函數在 x86-64 沒有 FMA 變體（遮掉 FMA 時 `cbrt` 0 個不同），推測是 aarch64 的 libm 在編譯時被自動 FMA 融合（例一的機制），而 x86-64 基本版是分開的乘、加。

### 例四：**同一個 binary、同一個 glibc，只是 CPU 不同**

glibc 在執行期依 CPU 功能以 IFUNC 選擇 `sin`/`cos`/`exp`/`log`/`pow`/`atan2`… 的 FMA 版本（Q2.2）。在同一台機器上用 `GLIBC_TUNABLES=glibc.cpu.hwcaps=-AVX2,-FMA` 模擬「沒有 FMA 的舊 CPU」：`sin` 有 686 個、`exp` 611、`pow` 726 個（約 0.07%）結果不同（1–2 ulp）。也就是說：**同一個執行檔，在 Haswell 之前與之後的 Intel CPU 上會算出不同的 `sin`**。

### 例五：核心數不同 → 平行歸約的結果不同

`#pragma omp parallel for reduction(+:s)` 的預設執行緒數 = CPU 核心數，每個執行緒加總自己的區段後再合併，捨入順序取決於執行緒數：

![](shots/q5_threads.png){.shot}

1、2、3、4 個執行緒得到 4 個不同的總和（差異 ~10⁻⁷，相對 10⁻¹³）。程式完全合法，但在 8 核和 64 核的機器上輸出不同。更進一步：4 執行緒的結果在兩次執行間也不同（本報告開發時一次得到 …164913、`run_all.sh` 中得到 …164908）—— libgomp 以執行緒**完成的先後順序**合併部分和，所以連同一台機器都不保證可重現。

**如何得到可重現的結果**：`-ffp-contract=off`（或只在需要時使用 `std::fma`）、不用 `long double`、使用正確捨入的數學函式庫（如 CORE-MATH）或自帶實作、固定歸約順序（固定分塊數而非依執行緒數分塊，或用 Kahan/exact accumulator）、在測試中以誤差容忍度而不是逐位元相等比較結果。

<div class="pagebreak"></div>

## Q6. 總結

### Q6.1 我對浮點數運算的看法

浮點數不是「不安全」，而是一個**規格嚴謹的近似系統**：IEEE 754 保證基本運算正確捨入，所以單一運算的誤差 ≤ ½ ulp，完全可預測。問題出在我們把它當成實數：

1. **誤差會累積與放大**：抵銷（相近數相減）、長鏈加總、病態問題都會放大捨入誤差 —— 需要分析或實驗量化，而不是祈禱。
2. **「同樣的程式」不保證同樣的位元**：編譯器（contraction、`-ffast-math`、向量化）、函式庫（libm 精度是 implementation-defined）、硬體（FMA dispatch、`long double` 格式）、平行化（歸約順序）都會改變結果（Q3、Q5）。
3. **通常截斷誤差比捨入誤差大**：本引擎中 RK4 的截斷誤差 ~10⁻⁸–10⁻⁹，double 的捨入誤差 ~10⁻¹⁵；RAPTOR 預設 RK2 的誤差甚至到 10⁻³。**先確保演算法收斂，再談型別。**
4. **驗證手段**：與解析解比較（陰影邊界）、收斂階數檢查（h⁴）、守恆量監控（E、L、Q、H）、以更高精度重跑（`__float128` 參考解）、與獨立實作（RAPTOR）交叉比對。

**實務準則**：預設 `double`；float 只在可向量化/受頻寬限制時使用；不用 `==` 比較浮點數（用相對+絕對容忍度）；不全域開 `-ffast-math`；小心 NaN/Inf 的傳播並在邊界檢查；輸出十進位時知道 printf 在做什麼；需要可重現時固定 contraction、libm 與歸約順序。

### Q6.2 額外細節：問題本身的「條件數」—— 再多位數也救不了病態問題

前面的題目都在談「計算」的誤差，但有一種誤差來源與浮點數格式無關：**問題本身對輸入的敏感度（condition number）**。黑洞陰影邊緣附近的光線會繞光子環好幾圈才離開，光子軌道是**不穩定**的：每繞半圈，擾動被放大 e<sup>γ</sup> 倍（γ 為 Lyapunov 指數）。`experiments/q6_extra.cpp` 讓兩條光線的撞擊參數只差 10⁻¹³（相對），量最終方向的差：

![](shots/q6.png){.shot}

![](fig/chaos.png){style="max-width:60%"}
<div class="caption">圖 10：誤差放大倍率 ∝ 1/δ（δ = 與陰影邊緣的相對距離）。δ = 10⁻¹⁰ 時放大 7×10⁹ 倍。</div>

在 δ = 10⁻¹⁰ 時，輸入 1 ulp 等級的差異（10⁻¹⁶）會變成 10⁻⁶ rad 的方向差；δ 再小，double 的捨入誤差就足以讓「捕獲或逃逸」的結論翻轉。**這不是 bug，也不能靠換 `long double` 解決**（只能把牆往後推幾位數），而是物理問題本身的性質 —— 這也解釋了第 1.4 節中 RAPTOR 與本引擎的最大誤差為什麼都集中在陰影邊緣。正確的做法是：辨識病態區域、在那裡以統計量（例如像素平均通量）而非單一光線來下結論。

同一個程式中另外兩個小陷阱（同樣在 `results/q6_extra.txt`）：

* **代數等價 ≠ 數值等價（catastrophic cancellation）**：RAPTOR 設定光線時計算 k<sub>θ</sub> = √(Q − L² cot² i + cos² i)，代入後代數上恆等於 β，但當 \|α\| ≫ \|β\| 時是大數相減：float 下 α = 20、β = 0.003 時算出 0.00478，**相對誤差 59%**。直接寫 k<sub>θ</sub> = β 就沒有這個問題（引擎提供 `exact_ktheta` 選項）。
* **非正規數很慢**：相依乘法鏈中運算元為非正規數時，每次乘法 45.7 ns，正常值 2.8 ns，**慢 17 倍**（microcode assist）—— 這也是 `-ffast-math` 設定 FTZ/DAZ 的動機，但代價是 Q3 的 #6。

<div class="pagebreak"></div>

## Q7. AI 使用自我揭露

> **（請依實際情況修改本段。）**

本作業使用 Anthropic 的 AI 程式助理 Claude（Claude Code）協助。

* **AI 完成的部分**：依我的指示（「寫一個 Kerr 度規測地線模擬引擎並與 RAPTOR 比較、加入平行化技術與快取命中率的比較」）撰寫了引擎程式碼、RAPTOR 橋接與比對程式、所有實驗程式、執行腳本、圖表、截圖工具與本報告初稿；並查閱了 glibc 2.39 的 `s_sin.c` 原始碼作為 Q2 的依據。
* **我自己做的部分**：（請填寫，例如：決定報告主題與實驗方向、閱讀並理解每一個實驗與程式碼、在自己的電腦上重跑 `run_all.sh` 核對數據、修改與補充說明、最終審閱。）
* 所有數據都來自實際執行的輸出（`results/*.txt`），截圖是把這些終端輸出以終端機樣式渲染成圖片（`scripts/screenshots.js`）。
