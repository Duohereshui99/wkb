# Stokes图分析：库仑+核势的共振结构

## Stokes Graph Analysis for Coulomb + Nuclear Potential

---

## 1. Exact WKB 基础

### 1.1 WKB解的形式

薛定谔方程：
```
d²ψ/dx² + Q(x)ψ = 0
```

其中：
```
Q(x) = (2μ/ℏ²)[E - V(x)] = k²(x)
```

**Exact WKB解**：
```
ψ_±(x) = 1/√S_odd(x) · exp(±∫ S_odd(x')dx')
```

其中 S_odd 是奇函数展开：
```
S_odd = S₀ + ℏ²S₂ + ℏ⁴S₄ + ...

S₀ = √Q(x) = k(x)
S₂ = -Q''/8Q^{3/2} + 5(Q')²/32Q^{5/2}
...
```

### 1.2 Stokes线定义

**Stokes线**是满足以下条件的复平面曲线：
```
Im ∫_{x_t}^{x} S_odd(x')dx' = 0
```

其中 x_t 是转折点（Q(x_t) = 0）。

**物理意义**：
- 在Stokes线上，ψ₊ 和 ψ₋ 的相对大小发生剧烈变化
- 穿越Stokes线时，次主导项可能变成主导项

### 1.3 Anti-Stokes线

**Anti-Stokes线**满足：
```
Re ∫_{x_t}^{x} S_odd(x')dx' = 0
```

**物理意义**：
- 在Anti-Stokes线上，|ψ₊| = |ψ₋|
- 两个解的幅度相等

---

## 2. 纯库仑势的Stokes图

### 2.1 库仑势的Q(r)

对于库仑势 V(r) = Z₁Z₂e²/r：
```
Q(r) = k² - 2ηk/r - (ℓ+1/2)²/r²
```

**转折点**：Q(r_t) = 0
```
k²r² - 2ηkr - (ℓ+1/2)² = 0

r_± = [η ± √(η² + (ℓ+1/2)²)] / k
```

### 2.2 Stokes图结构（s波，ℓ=0）

对于 E > 0（散射态）：

```
                    Im r
                      │
                      │     Stokes线 (指向 Im r → +∞)
                      │    ╱
                      │   ╱
                      │  ╱
                      │ ╱
    ──────────────────┼─────●──────────────→ Re r
                      │    r₊
                      │ ╲
                      │  ╲
                      │   ╲
                      │    ╲
                      │     Stokes线 (指向 Im r → -∞)
```

**关键特征**：
- 只有一个物理转折点 r₊ > 0
- Stokes线从 r₊ 出发，延伸到 |Im r| → ∞
- 没有经典禁戒区（无隧穿）

### 2.3 Stokes图结构（高ℓ或低能）

当 η 较大或 ℓ > 0 时，可能出现两个转折点：

```
                    Im r
                      │
                (-)   │   (+)
                      │ ╲ ╱
    ──────────────────┼──╳────●────────→ Re r
                      │ ╱r₁╲ r₂
                (+)   │     (-)
                      │
```

**经典禁戒区**：r₁ < r < r₂（势垒区域）

---

## 3. 库仑+核势的Stokes图

### 3.1 势能结构

```
V_eff(r) = V_N(r) + V_C(r) + V_cent(r)
```

典型形状：

```
V(r)
  │
  │         ╭──╮  Coulomb barrier
  │        ╱    ╲
E ├───────╱──────╲─────────────────
  │   ┌──╯        ╲
  │   │            ╲______________
  │   │   Nuclear
  │   │   pocket
  │   └──────┐
  │          V_N
  └────────────────────────────────→ r
      r₁   r₂        r₃
```

### 3.2 三转折点情况

对于 0 < E < V_barrier：

**转折点**：
- r₁：核势阱的外边缘
- r₂：库仑势垒的内侧
- r₃：库仑势垒的外侧

**区域划分**：
```
区域 I:   r < r₁     (核内，振荡)      Q > 0
区域 II:  r₁ < r < r₂ (势阱-势垒过渡)  Q < 0 或 Q > 0
区域 III: r₂ < r < r₃ (势垒内，指数)    Q < 0
区域 IV:  r > r₃     (渐近区，振荡)    Q > 0
```

### 3.3 复平面Stokes图

**完整的Stokes图**（能量低于势垒）：

```
                         Im r
                           │
                           │
            Stokes线       │      Stokes线
           (从r₁出发)     │     (从r₃出发)
              ╲           │         ╱
               ╲          │        ╱
                ╲         │       ╱
    ─────────────●────────┼───●──●─────────→ Re r
                r₁        │  r₂ r₃
                 ╲        │       ╲
                  ╲       │        ╲
                   ╲      │         ╲
                          │
                          │
```

**共振能量时的Stokes图变化**：

当能量接近准束缚态能量时，Stokes图发生拓扑变化：
- 某些Stokes线可能连接不同的转折点
- 形成新的拓扑结构，对应共振态

---

## 4. Siegert边界条件的几何意义

### 4.1 散射态 vs 共振态

**散射态**：
- 能量为实数 E ∈ ℝ
- 波函数在 r → ∞ 振荡
- 同时存在入射波和出射波

**共振态**：
- 能量为复数 E = E_R - iΓ/2
- 波函数在实轴上发散
- 只有出射波（Siegert条件）

### 4.2 Siegert条件的Stokes图表述

**Siegert边界条件**：在 r → ∞ 只有出射波

在Exact WKB框架中：
```
ψ⁻ 在 Im r → -∞ 处占主导（对于 Re r > 0）
```

**几何意义**：
```
                    Im r
                      │
                      │
                      │ 出射波主导
    ──────────────────┼──────●────────→ Re r
                      │     r₃
                      │╲
                      │ ╲
                      │  ╲  Siegert路径
                      │   ╲ (沿Stokes线延伸)
                      │    ╲
                      │     ● 共振波函数可归一化
```

### 4.3 与S矩阵极点的关系

**散射矩阵**：
```
S(k) = e^{2iδ(k)}
```

**极点条件**：
```
1/S(k_R) = 0  ⟺  e^{-2iδ(k_R)} = 0  ⟺  δ(k_R) → +i∞
```

**物理解释**：
- 入射波振幅为零
- 只有出射波
- 对应Siegert边界条件

---

## 5. 共振能量的量子化

### 5.1 Bohr-Sommerfeld条件

**准束缚态能量**由以下条件确定：
```
∮ k(r) dr = 2πℏ(n + 1/2)
```

对于核势阱中的准束缚态：
```
∫_{r₁}^{r₂} √[2μ(E - V_eff)/ℏ²] dr = (n + 1/2)π
```

### 5.2 Exact WKB量子化

**精确量子化条件**（包含所有ℏ修正）：
```
∮_C S_odd(r) dr = 2πiℏ(n + 1/2)
```

其中 C 是围绕两个转折点的闭合路径。

**Borel重求和**：
```
S_odd^{Borel}(r) = ∫_0^∞ e^{-t} Ŝ_odd(r, ℏt) dt
```

处理发散的WKB级数。

### 5.3 共振宽度

**Gamow公式**：
```
Γ = (ℏω/2π) · T
```

其中：
- ω = 核内振动频率 ≈ v/(2R)
- T = exp(-2γ) 是穿透因子
- γ = (1/ℏ)∫_{r₂}^{r₃} |p(r)| dr 是Gamow因子

**Exact WKB修正**：
```
Γ = Γ_WKB · [1 + O(ℏ²)]
```

包含高阶量子修正。

---

## 6. 与库仑波函数的连接

### 6.1 渐近匹配

在 r → ∞，波函数渐近为库仑波函数：
```
ψ(r → ∞) → A·F_ℓ(η,kr) + B·G_ℓ(η,kr)
```

**WKB表示**：
```
ψ_WKB(r → ∞) → 1/√k · sin(kr - ηln(2kr) - ℓπ/2 + σ_ℓ + δ_ℓ)
```

**匹配条件**：
```
tan(δ_ℓ) = -B/A
```

### 6.2 S矩阵元

**库仑修正的S矩阵**：
```
S_ℓ = e^{2i(σ_ℓ + δ_ℓ)} = e^{2iσ_ℓ} · (A + iB)/(A - iB)
```

**WKB近似**：
```
δ_ℓ^{WKB} = lim_{r→∞}[∫_{r_t}^{r} k_{eff}(r')dr' - kr + ηln(2kr)] - (ℓ+1/2)π/2
```

### 6.3 极点位置

**共振位置**由以下方程确定：
```
A(k_R) - iB(k_R) = 0
```

或等价地：
```
cot(δ_ℓ(k_R)) = i
```

在复k平面：
```
k_R = k_r - i·k_i,  其中 k_r > 0, k_i > 0
```

对应复能量：
```
E_R = ℏ²k_R²/(2μ) = (ℏ²k_r²/2μ - ℏ²k_i²/2μ) - i(ℏ²k_r·k_i/μ)
     = E_res - iΓ/2
```

---

## 7. 数值方法

### 7.1 Stokes图的数值绘制

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def plot_stokes_graph(V_eff, E, mu, hbar, r_range, n_points=100):
    """
    绘制Stokes图

    参数:
        V_eff: 有效势函数 V(r)
        E: 能量
        mu: 约化质量
        hbar: 约化普朗克常数
        r_range: (r_min, r_max, im_max)
        n_points: 网格点数
    """
    r_min, r_max, im_max = r_range

    # 创建复平面网格
    re_r = np.linspace(r_min, r_max, n_points)
    im_r = np.linspace(-im_max, im_max, n_points)
    RE, IM = np.meshgrid(re_r, im_r)
    R_complex = RE + 1j * IM

    # 计算 Q(r) = k²(r)
    def Q(r):
        if np.real(r) <= 0:
            return np.inf
        return 2 * mu * (E - V_eff(np.real(r))) / hbar**2

    # 找转折点
    def find_turning_points():
        r_test = np.linspace(r_min, r_max, 1000)
        Q_vals = np.array([Q(r) for r in r_test])
        # 找符号变化点
        sign_changes = np.where(np.diff(np.sign(Q_vals)))[0]
        return r_test[sign_changes]

    turning_points = find_turning_points()

    # 对每个转折点，计算从该点出发的积分
    fig, ax = plt.subplots(figsize=(10, 8))

    for r_t in turning_points:
        # 计算相位积分 ∫ S_odd dr
        phase_real = np.zeros_like(RE)
        phase_imag = np.zeros_like(RE)

        for i in range(n_points):
            for j in range(n_points):
                r = R_complex[i, j]
                # 简化：只用S_0 = √Q
                # 实际应该用完整的S_odd
                try:
                    # 沿实轴积分的近似
                    re_part = np.real(r)
                    if re_part > r_t:
                        integral, _ = quad(lambda x: np.sqrt(abs(Q(x))),
                                          r_t, re_part)
                        phase_real[i, j] = integral
                        phase_imag[i, j] = np.imag(r) * np.sqrt(abs(Q(re_part)))
                except:
                    pass

        # 绘制Stokes线 (Im∫S_odd = 0)
        ax.contour(RE, IM, phase_imag, levels=[0], colors='blue',
                  linewidths=2, linestyles='-')

        # 绘制Anti-Stokes线 (Re∫S_odd = 0)
        ax.contour(RE, IM, phase_real, levels=[0], colors='red',
                  linewidths=1, linestyles='--')

        # 标记转折点
        ax.plot(r_t, 0, 'ko', markersize=10)
        ax.annotate(f'$r_t={r_t:.1f}$', (r_t, 0.5))

    ax.set_xlabel('Re(r)', fontsize=12)
    ax.set_ylabel('Im(r)', fontsize=12)
    ax.set_title('Stokes Graph', fontsize=14)
    ax.axhline(y=0, color='gray', linewidth=0.5)
    ax.axvline(x=0, color='gray', linewidth=0.5)
    ax.grid(True, alpha=0.3)

    # 图例
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='blue', linewidth=2, label='Stokes line'),
        Line2D([0], [0], color='red', linewidth=1, linestyle='--',
               label='Anti-Stokes line'),
        Line2D([0], [0], marker='o', color='black', linestyle='',
               markersize=10, label='Turning point')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    return fig
```

### 7.2 共振搜索算法

```python
def find_resonance_exact_wkb(wkb_system, E_range, ell=0):
    """
    用Exact WKB方法搜索共振

    算法:
    1. 在实能量轴上计算相移
    2. 找 δ = π/2 的能量（共振中心）
    3. 用Gamow公式计算宽度
    4. 或直接在复能量平面搜索极点
    """
    from scipy.optimize import brentq

    # 步骤1: 扫描相移
    E_vals = np.linspace(E_range[0], E_range[1], 100)
    delta_vals = np.array([wkb_system.wkb_phase_shift(E, ell) for E in E_vals])

    # 步骤2: 找共振能量 (δ = π/2 + nπ)
    resonances = []

    for n in range(5):  # 搜索前几个共振
        target = np.pi/2 + n * np.pi

        # 找穿越点
        for i in range(len(E_vals) - 1):
            if (delta_vals[i] - target) * (delta_vals[i+1] - target) < 0:
                # 精确定位
                E_res = brentq(
                    lambda E: wkb_system.wkb_phase_shift(E, ell) - target,
                    E_vals[i], E_vals[i+1]
                )

                # 步骤3: 计算宽度
                gamma, T = wkb_system.gamow_factor(E_res, ell)
                Gamma = wkb_system.decay_width(E_res, ell)[0]

                resonances.append({
                    'n': n,
                    'E_res': E_res,
                    'Gamma': Gamma,
                    'gamow_factor': gamma
                })
                break

    return resonances
```

---

## 8. 总结

### 8.1 方法对比

| 方法 | 优点 | 缺点 |
|-----|------|------|
| 数值S矩阵极点 | 精确 | 计算量大，需要复能量延拓 |
| 标准WKB | 简单直观 | 低阶近似，转折点处不准确 |
| Exact WKB | 系统精确 | 理论复杂，需要Borel重求和 |
| ERF方法 | 解析延拓方便 | 依赖低能展开参数 |

### 8.2 推荐工作流程

1. **初步估计**：用标准WKB/Gamow公式
2. **精确计算**：用数值S矩阵极点方法
3. **解析理解**：用Exact WKB分析Stokes图结构
4. **参数提取**：用ERF方法得到散射长度等可观测量

### 8.3 未来方向

- 将Exact WKB推广到包含库仑长程势
- 发展库仑修正的连接公式
- 与ab initio核结构计算结合

---

*文档版本: 1.0*
*最后更新: 2025-12-22*
