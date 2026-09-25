# NTHU-AdvCPP-HW2
# Floating Point 行為專案報告

### Advanced C++: Large-Scale Projects and Object-Oriented Design Homework 2

**Deadline** : 2026/10/07 23:59 前於 eeclass 上繳交 PDF 格式的報告

---

> **報告情境：**
> 
> 假設現在你正打算撰寫一個電腦程式用於科學計算、或是開發遊戲需要處裡物理碰撞等等，導致你想使用浮點數來計算一些數值。然而你的主管或是專案的其他隊友，知道浮點數計算是很不安全的。你必須說服這些人，你知道浮點數的計算怎樣是安全的，怎樣是危險的。
> 
> 你需要設計一些實驗，來評估：
> 1. 對於相同的數值，在不同浮點數型態 (float / double …) 計算上的效能、精準度等差異，提出你的 Benchmark
> 2. 詳細的介紹浮點數的硬體實作
> 
> 報告請按照助教提出的問題，作為整個文件的大綱來完成。並對關鍵結果，給出你的執行輸出截圖來佐證。
> 
> **繳交格式要求：**
> 
> 繳交的報告，需要命名為 : **學號_姓名.pdf** (如 `112000000_王大明.pdf`)。並且第一頁，需要有下列資訊作為標題：
> 
> <div align="center">
> 
> Advanced C++ Homework 2  
> 學號, 姓名  
> 
> </div>

---

### Q0. 實驗環境
你可以選擇任意的一個能進行浮點數計算的系統來進行你的實驗，不限於 PC、筆電、手機、遊戲機。  
但是你需要能清楚地描述實驗的規格 (什麼硬體、什麼編譯器等等)，達到「能夠被閱讀報告的第三方重現你實驗」的程度。

### Q1. Benchmark
請提出對於 C++ 可用的小數型態 (float、double、long double、或是 float128 等等) 的基本比較，計算速度的差異。  
不僅是列出實驗數據，進一步的給出為什麼有差異的解釋。

### Q2. Library Implementation
詳細的介紹 `cmath` 中，自選**一個**函數 (如：`sin(x)`) ，在 C++ 上的規格是什麼，與具體如何被實作，實作部分請指出其關鍵演算法與實作細節即可，例如：如何縮小輸入範圍、如何進行數值近似、是否使用特殊的硬體指令等。不需要逐行介紹函式庫的原始碼。

### Q3. IEEE 754 與 -ffast-math
GCC 編譯器存在一個叫做 `-ffast-math` 的參數 ，這有什麼用？詳細的介紹關閉與開啟此選項的差異，並給出與解釋一個僅因為有無開啟 `-ffast-math` 選項 ，而導致計算結果發生改變的例子。

### Q4. 浮點數輸出
```cpp
printf("%.9f") 
std::fixed << setprecision(9)
```
這樣寫真的是四捨五入 (到小數後第 9 位) 嗎？請判斷這個敘述是否精確，並提出反例。

### Q5. 加分題
找出一個滿足以：**相同的程式碼、相同版本的編譯器、相同的編譯參數，沒有未定義行為，**但是浮點數計算結果可能因為硬體、作業系統、ABI 或其他執行環境的差異而導致結果不同的例子。  
本題目可以與其他同學合作。

### Q6. 總結
給出你對浮點數運算的看法，有什麼需要注意的地方，總結整份報告，並提出一個在上面問題中未提到，但是也是浮點數計算需要注意的細節。

### Q7. AI 使用自我揭露 (不記入作業分數)
簡要描述你自己做了什麼，請 AI 幫你做了什麼。

---

## 實作：Kerr 黑洞測地線模擬引擎 + 浮點數實驗

報告以「黑洞光線追跡 (general-relativistic ray tracing)」這個真實的科學計算問題作為主軸：
自行實作一個 C++ Kerr 時空零測地線 (null geodesic) 模擬引擎，並與天文界實際使用的
[RAPTOR](https://github.com/tbronzwaer/raptor)（Bronzwaer et al. 2018, A&A 613, A2）逐條光線比對，
再用這個引擎與一系列小實驗回答 Q0–Q7。報告本體：[`report/report.md`](report/report.md)（PDF 由 `scripts/build_report.py` 產生）。

### 目錄

| 路徑 | 內容 |
|---|---|
| `include/kerr/` | header-only 引擎：`real_math.hpp`（float/double/long double/`__float128`/`std::float16_t`/SIMD 的數學函式包裝）、`kerr.hpp`（Kerr–Schild 座標下的 Hamiltonian 測地線方程、相機、薄盤）、`tracer.hpp`（RK2/RK4/Dormand–Prince 5(4) 積分器）、`render.hpp`、`critical_curve.hpp`（Bardeen 解析陰影邊界） |
| `src/kerr_render.cpp` | 產生黑洞影像（可選浮點型別、積分器） |
| `src/kerr_validate.cpp` | 與解析解（Kerr 陰影邊界）比對，驗證收斂階數與各型別精度極限 |
| `src/raptor_compare.cpp`, `raptor_bridge/` | 以 RAPTOR 原始碼（未修改）逐條光線與本引擎比對 |
| `src/kerr_precision.cpp` | Q1：同一個光線追跡工作量在各浮點型別下的速度與 round-off |
| `src/parallel_bench.cpp` | 平行化技術比較：std::thread / OpenMP 各種 schedule / `std::execution::par` / `std::experimental::simd` 光線封包 / false sharing |
| `experiments/` | Q1 micro-benchmark、快取實驗、Q2 `sin`、Q3 `-ffast-math`、Q4 `printf`、Q5 跨平台、Q6 額外細節 |
| `scripts/` | `fetch_raptor.sh`、`run_all.sh`（跑全部實驗 → `results/`）、`make_figures.py`、`screenshots.js`、`build_report.py` |
| `results/` | 所有實驗的原始終端輸出（報告中的數據都出自這裡） |

### 重現

```bash
# Ubuntu 24.04: build-essential libgsl-dev libmpfr-dev libgmp-dev libtbb-dev valgrind
#               g++-aarch64-linux-gnu qemu-user musl-tools python3-numpy python3-matplotlib
make raptor              # 下載 RAPTOR（固定 commit 08cb9a2）
make all cross           # 編譯全部（含 aarch64 / musl 版本）
./scripts/run_all.sh     # 約 30 分鐘，輸出到 results/
python3 scripts/make_figures.py && node scripts/screenshots.js && python3 scripts/build_report.py
```
