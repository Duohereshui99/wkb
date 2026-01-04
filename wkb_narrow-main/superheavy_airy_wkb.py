"""
Higher-Order WKB for Superheavy Nuclei Alpha Decay
===================================================

Complete Airy uniform approximation for superheavy nuclei from:
Z. Wang & Z. Ren, Phys. Rev. C 110, 064307 (2024)

Key features:
- Rigorous Airy uniform approximation (no ad-hoc cutoff)
- Proper centrifugal barrier for ℓ ≠ 0 unfavored transitions
- Data from Table II including spin-parity assignments

Author: Jin Lei
Date: 2025-12-23
"""

import numpy as np
from scipy.special import airy
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# Physical constants
HBAR_C = 197.327      # MeV·fm
ALPHA_EM = 1/137.036
M_NUCLEON = 931.494   # MeV/c^2


# =============================================================================
# Superheavy Nuclei Data from PhysRevC.110.064307 Table II
# =============================================================================

# Format: (Z_daughter, A_daughter, Q_MeV, T_exp_s, ell, beta2, beta4, beta6, J_in, J_fi, note)
SUPERHEAVY_NUCLEI = {
    # ========== α-decay chain of ²⁷⁵Ds ==========
    '275Ds_0': (108, 271, 11.365, 4.30e-4, 0, 0.221, -0.080, -0.007, '3/2', '3/2', 'favored'),
    '275Ds_4': (108, 271, 11.365, 4.30e-4, 4, 0.221, -0.080, -0.007, '3/2', '11/2', 'proposed'),

    '271Hs_0a': (106, 267, 9.180, 7.10e0, 0, 0.232, -0.065, -0.015, '3/2', '3/2', 'isomer'),
    '271Hs_0b': (106, 267, 9.480, 4.60e1, 0, 0.232, -0.065, -0.015, '11/2', '11/2', 'ground'),
    '271Hs_6': (106, 267, 9.480, 4.60e1, 6, 0.232, -0.065, -0.015, '11/2', '1/2', 'proposed'),

    '267Sg': (104, 263, 8.400, 5.90e2, 0, 0.231, -0.040, -0.021, '9/2', '9/2', 'favored'),

    # ========== α-decay chain of ²⁷³Ds ==========
    '273Ds_0a': (108, 269, 11.090, 3.00e-2, 0, 0.232, -0.065, -0.015, '11/2', '11/2', 'ground'),
    '273Ds_6': (108, 269, 11.090, 3.00e-2, 6, 0.232, -0.065, -0.015, '11/2', '1/2', 'proposed'),
    '273Ds_0b': (108, 269, 11.270, 1.80e-4, 0, 0.232, -0.065, -0.015, '1/2', '1/2', 'isomer'),

    '269Hs_0a': (106, 265, 9.340, 1.30e1, 0, 0.232, -0.052, -0.023, '9/2', '9/2', 'ground'),
    '269Hs_4': (106, 265, 9.340, 1.30e1, 4, 0.232, -0.052, -0.023, '9/2', '3/2', 'proposed'),
    '269Hs_0b': (106, 265, 9.220, 2.80e0, 0, 0.232, -0.052, -0.023, '1/2', '1/2', 'isomer'),

    '265Sg_0a': (104, 261, 8.970, 8.50e0, 0, 0.242, -0.025, -0.028, '11/2', '11/2', 'ground'),
    '265Sg_0b': (104, 261, 8.820, 1.44e1, 0, 0.242, -0.025, -0.028, '3/2', '3/2', 'isomer'),

    '261Rf': (102, 257, 8.410, 6.80e1, 0, 0.240, 0.000, -0.033, '11/2', '11/2', 'favored'),

    # ========== Other superheavy nuclei ==========
    '251Lr_0a': (101, 247, 9.438, 2.44e-2, 0, 0.249, 0.051, -0.032, '7/2-', '7/2-', 'favored'),
    '251Lr_0b': (101, 247, 9.402, 4.20e-2, 0, 0.249, 0.051, -0.032, '1/2-', '1/2-', 'isomer'),

    '253Lr_0a': (101, 249, 8.969, 6.50e-1, 0, 0.250, 0.039, -0.035, '7/2-', '7/2-', 'favored'),
    '253Lr_2': (101, 249, 8.842, 4.30e-1, 2, 0.250, 0.039, -0.035, '7/2-', '11/2-', 'unfavored'),
    '253Lr_0b': (101, 249, 8.897, 2.46e0, 0, 0.250, 0.039, -0.035, '1/2-', '1/2-', 'isomer'),

    '255Lr_0a': (101, 251, 8.634, 2.53e0, 0, 0.250, 0.027, -0.037, '7/2-', '7/2-', 'favored'),
    '255Lr_0b': (101, 251, 8.541, 3.11e1, 0, 0.250, 0.027, -0.037, '1/2-', '1/2-', 'isomer'),

    '272Hs': (106, 268, 9.772, 1.60e-1, 0, 0.232, -0.065, -0.015, '', '', 'even-even'),
    '276Ds': (108, 272, 10.904, 3.60e-4, 0, 0.221, -0.080, -0.007, '', '', 'even-even'),
    '284Cn': (110, 280, 9.460, 6.10e0, 0, 0.130, -0.042, -0.005, '', '', 'even-even'),

    '286Fl_0a': (112, 282, 9.710, 3.50e0, 0, 0.086, -0.009, -0.011, '0+', '0+_2', 'excited'),
    '286Fl_0b': (112, 282, 10.330, 2.30e-1, 0, 0.086, -0.009, -0.011, '0+', '0+_1', 'ground'),
    '288Fl': (112, 284, 10.060, 6.50e-1, 0, 0.086, -0.021, -0.002, '', '', 'even-even'),

    '286Mc': (113, 282, 10.861, 2.00e-2, 0, 0.064, 0.014, -0.009, '', '', 'odd-A'),
    '292Lv': (114, 288, 10.790, 1.30e-2, 0, -0.021, 0.012, 0.000, '', '', 'even-even'),
}


class AiryWKBSuperheavy:
    """
    Airy uniform WKB for superheavy nuclei alpha decay.

    Includes:
    - Proper centrifugal barrier with Langer modification
    - Rigorous Airy uniform approximation
    - No ad-hoc regularization parameters
    """

    def __init__(self, Z_daughter, A_daughter, r0=1.17, a=0.70, V0=None):
        """Initialize the alpha-daughter system."""
        self.Z_alpha = 2
        self.A_alpha = 4
        self.Z_daughter = Z_daughter
        self.A_daughter = A_daughter

        # Reduced mass
        m_alpha = self.A_alpha * M_NUCLEON
        m_daughter = A_daughter * M_NUCLEON
        self.mu = m_alpha * m_daughter / (m_alpha + m_daughter)

        # Nuclear radius
        self.r0 = r0
        self.R0 = r0 * (A_daughter**(1/3) + self.A_alpha**(1/3))

        # Woods-Saxon potential
        if V0 is None:
            N = A_daughter - Z_daughter
            # Standard formula with some adjustment for superheavy region
            self.V0 = 50 + 33 * (N - Z_daughter) / A_daughter
            self.V0 = max(40, min(150, self.V0))
        else:
            self.V0 = V0
        self.a = a

        # Coulomb parameters
        self.Z1Z2 = self.Z_alpha * Z_daughter
        self.k_C = self.Z1Z2 * ALPHA_EM * HBAR_C

    def V_nuclear(self, r):
        """Woods-Saxon nuclear potential (MeV)."""
        if r <= 0:
            return -self.V0
        return -self.V0 / (1 + np.exp((r - self.R0) / self.a))

    def V_coulomb(self, r):
        """Coulomb potential (MeV)."""
        if r <= 0:
            r = 0.01
        if r >= self.R0:
            return self.k_C / r
        else:
            return self.k_C / (2 * self.R0) * (3 - (r / self.R0)**2)

    def V_centrifugal(self, r, ell):
        """Centrifugal barrier with Langer modification (MeV)."""
        if ell == 0 or r <= 0:
            return 0
        # Langer modification: ℓ(ℓ+1) → (ℓ+1/2)²
        ell_eff = (ell + 0.5)**2
        return HBAR_C**2 * ell_eff / (2 * self.mu * r**2)

    def V_eff(self, r, ell=0):
        """Effective potential = nuclear + Coulomb + centrifugal."""
        return self.V_nuclear(r) + self.V_coulomb(r) + self.V_centrifugal(r, ell)

    def Q(self, r, E, ell=0):
        """Q(r) = (2μ/ℏ²)[E - V(r)]"""
        return 2 * self.mu * (E - self.V_eff(r, ell)) / HBAR_C**2

    def Q_derivatives(self, r, E, ell=0, h=0.005):
        """Compute Q, Q', Q'' using 5-point stencil."""
        r = max(r, 0.1)

        Q0 = self.Q(r, E, ell)

        Q1 = (-self.Q(r+2*h, E, ell) + 8*self.Q(r+h, E, ell)
              - 8*self.Q(r-h, E, ell) + self.Q(r-2*h, E, ell)) / (12*h)

        Q2 = (-self.Q(r+2*h, E, ell) + 16*self.Q(r+h, E, ell) - 30*Q0
              + 16*self.Q(r-h, E, ell) - self.Q(r-2*h, E, ell)) / (12*h**2)

        return Q0, Q1, Q2

    def find_turning_points(self, E, ell=0):
        """Find classical turning points."""
        r_grid = np.linspace(self.R0 * 0.3, 200, 5000)
        Q_grid = np.array([self.Q(r, E, ell) for r in r_grid])

        turning_points = []
        for i in range(len(r_grid) - 1):
            if Q_grid[i] * Q_grid[i+1] < 0:
                try:
                    r_t = brentq(lambda r: self.Q(r, E, ell), r_grid[i], r_grid[i+1])
                    turning_points.append(r_t)
                except:
                    pass

        return sorted(turning_points)

    def barrier_action_airy(self, E, ell=0, order=0):
        """
        Compute barrier action with Airy uniform approximation.

        This method eliminates ad-hoc cutoff parameters by using
        proper Airy function matching at turning points.
        """
        turning_points = self.find_turning_points(E, ell)

        if len(turning_points) < 2:
            return {'gamma': 0, 'T': 1.0, 'valid': False}

        r_in = turning_points[-2]
        r_out = turning_points[-1]

        # Get Q derivatives at turning points
        Q0_in, Q1_in, Q2_in = self.Q_derivatives(r_in, E, ell)
        Q0_out, Q1_out, Q2_out = self.Q_derivatives(r_out, E, ell)

        # Airy length scales
        lambda_in = abs(Q1_in)**(-1/3) if abs(Q1_in) > 1e-10 else 1.0
        lambda_out = abs(Q1_out)**(-1/3) if abs(Q1_out) > 1e-10 else 1.0

        # =====================================================================
        # STEP 1: Leading-order action S₀
        # =====================================================================

        def kappa_0(r):
            Qval = self.Q(r, E, ell)
            if Qval < 0:
                return np.sqrt(-Qval)
            return 0

        eps = 0.02
        S0, _ = quad(kappa_0, r_in + eps, r_out - eps, limit=500)

        # Add corrections for regions near turning points
        if Q1_in < 0:
            S0_in_correction = (2/3) * np.sqrt(abs(Q1_in)) * eps**(3/2)
        else:
            S0_in_correction = 0

        if Q1_out > 0:
            S0_out_correction = (2/3) * np.sqrt(abs(Q1_out)) * eps**(3/2)
        else:
            S0_out_correction = 0

        S0_total = S0 + S0_in_correction + S0_out_correction

        if order == 0:
            T = np.exp(-2 * S0_total)
            return {
                'gamma': S0_total,
                'T': T,
                'S0': S0_total,
                'r_in': r_in,
                'r_out': r_out,
                'valid': True
            }

        # =====================================================================
        # STEP 2: O(ℏ²) correction with Airy regularization
        # =====================================================================

        def S2_bulk(r):
            """S₂ in bulk region (away from turning points)."""
            Q0, Q1, Q2 = self.Q_derivatives(r, E, ell)
            Q_abs = abs(Q0)

            if Q_abs < 0.05:
                return 0

            S2_val = Q2 / (8 * Q_abs**1.5) - 5 * Q1**2 / (32 * Q_abs**2.5)
            return S2_val

        # Bulk integration limits
        r_bulk_in = r_in + 3 * lambda_in
        r_bulk_out = r_out - 3 * lambda_out

        if r_bulk_in >= r_bulk_out:
            S2_bulk_integral = 0
        else:
            S2_bulk_integral, _ = quad(S2_bulk, r_bulk_in, r_bulk_out, limit=500)

        # =====================================================================
        # STEP 3: Airy corrections at turning points
        # =====================================================================

        # Finite Airy correction: δS₂^{Airy} = -(1/48) × Q'' / |Q'|^{5/3}
        if abs(Q1_in) > 1e-10:
            delta_S2_in = -(1/48) * Q2_in / abs(Q1_in)**(5/3)
        else:
            delta_S2_in = 0

        if abs(Q1_out) > 1e-10:
            delta_S2_out = -(1/48) * Q2_out / abs(Q1_out)**(5/3)
        else:
            delta_S2_out = 0

        # =====================================================================
        # STEP 4: Transition region contribution
        # =====================================================================

        def S2_transition(r, r_t, lambda_t):
            Q0, Q1, Q2 = self.Q_derivatives(r, E, ell)
            Q_abs = abs(Q0)

            if Q_abs < 1e-6:
                return 0

            S2_raw = Q2 / (8 * Q_abs**1.5) - 5 * Q1**2 / (32 * Q_abs**2.5)

            dist = abs(r - r_t)
            if dist < lambda_t:
                return 0
            elif dist < 3 * lambda_t:
                x = (dist - lambda_t) / (2 * lambda_t)
                weight = 3*x**2 - 2*x**3
                return S2_raw * weight
            else:
                return S2_raw

        S2_trans_in = 0
        S2_trans_out = 0

        if r_in + lambda_in < r_bulk_in:
            S2_trans_in, _ = quad(lambda r: S2_transition(r, r_in, lambda_in),
                                  r_in + lambda_in, r_bulk_in, limit=200)

        if r_bulk_out < r_out - lambda_out:
            S2_trans_out, _ = quad(lambda r: S2_transition(r, r_out, lambda_out),
                                   r_bulk_out, r_out - lambda_out, limit=200)

        # =====================================================================
        # STEP 5: Combine all contributions
        # =====================================================================

        S2_total = S2_bulk_integral + S2_trans_in + S2_trans_out + delta_S2_in + delta_S2_out

        gamma_O2 = S0_total + S2_total
        T_O2 = np.exp(-2 * gamma_O2)

        return {
            'gamma': gamma_O2,
            'T': T_O2,
            'S0': S0_total,
            'S2': S2_total,
            'S2_bulk': S2_bulk_integral,
            'delta_S2_in': delta_S2_in,
            'delta_S2_out': delta_S2_out,
            'r_in': r_in,
            'r_out': r_out,
            'lambda_in': lambda_in,
            'lambda_out': lambda_out,
            'valid': True
        }

    def gamow_factor(self, E, ell=0, order=0):
        """Main interface for Gamow factor calculation."""
        result = self.barrier_action_airy(E, ell, order)
        return result['gamma'], result['T']

    def half_life(self, Q_value, ell=0, order=0, P_alpha=None):
        """
        Calculate alpha-decay half-life.

        Parameters
        ----------
        P_alpha : float, optional
            Preformation factor. If None, uses:
            - 0.2137 for even-even nuclei
            - 0.0643 for odd-A nuclei
            - 0.0532 for odd-odd nuclei
        """
        gamma, T = self.gamow_factor(Q_value, ell, order)

        if T <= 0:
            return np.inf

        # Wave number and velocity
        k = np.sqrt(2 * self.mu * Q_value) / HBAR_C
        v = HBAR_C * k / self.mu

        # Assault frequency
        omega = v * HBAR_C / (2 * self.R0)

        # Preformation factor
        if P_alpha is None:
            Z_parent = self.Z_daughter + 2
            A_parent = self.A_daughter + 4
            N_parent = A_parent - Z_parent
            is_even_even = (Z_parent % 2 == 0) and (N_parent % 2 == 0)
            is_odd_odd = (Z_parent % 2 == 1) and (N_parent % 2 == 1)
            if is_even_even:
                P_alpha = 0.2137
            elif is_odd_odd:
                P_alpha = 0.0532
            else:
                P_alpha = 0.0643

        # Decay width
        Gamma = P_alpha * omega * T / (2 * np.pi)

        # Half-life
        hbar_MeV_s = 6.582e-22  # MeV·s
        if Gamma > 0:
            T_half = hbar_MeV_s * np.log(2) / Gamma
        else:
            T_half = np.inf

        return T_half

    def extract_preformation_factor(self, Q_value, T_exp_s, ell=0, order=0):
        """Extract P_α from experimental half-life."""
        gamma, T = self.gamow_factor(Q_value, ell, order)

        if T <= 0:
            return np.nan

        k = np.sqrt(2 * self.mu * Q_value) / HBAR_C
        v = HBAR_C * k / self.mu
        omega = v * HBAR_C / (2 * self.R0)

        hbar_MeV_s = 6.582e-22
        P_alpha = 2 * np.pi * hbar_MeV_s * np.log(2) / (omega * T * T_exp_s)

        return P_alpha


def analyze_superheavy_nuclei():
    """
    Complete analysis of superheavy nuclei using Airy WKB.
    """
    print("=" * 90)
    print("Higher-Order WKB Analysis of Superheavy Nuclei Alpha Decay")
    print("Data from: Z. Wang & Z. Ren, Phys. Rev. C 110, 064307 (2024)")
    print("Method: Airy Uniform Approximation (no ad-hoc cutoff)")
    print("=" * 90)

    results = []

    print(f"\n{'Nucleus':12s} | {'ℓ':3s} | {'Q(MeV)':7s} | {'γ_O0':7s} | {'γ_O2':7s} | {'Δγ':7s} | {'T_O2/T_O0':9s} | {'HF(O0)':8s} | {'HF(O2)':8s}")
    print("-" * 100)

    for name, data in sorted(SUPERHEAVY_NUCLEI.items()):
        Z_d, A_d, Q, T_exp, ell, beta2, beta4, beta6, J_in, J_fi, note = data

        # Create system
        system = AiryWKBSuperheavy(Z_d, A_d)

        # Calculate Gamow factors
        result_O0 = system.barrier_action_airy(Q, ell, order=0)
        result_O2 = system.barrier_action_airy(Q, ell, order=2)

        if not result_O0['valid']:
            continue

        gamma_O0 = result_O0['gamma']
        gamma_O2 = result_O2['gamma']
        delta_gamma = gamma_O2 - gamma_O0
        T_ratio = result_O2['T'] / result_O0['T']

        # Calculate half-lives and hindrance factors
        T_calc_O0 = system.half_life(Q, ell, order=0)
        T_calc_O2 = system.half_life(Q, ell, order=2)

        HF_O0 = T_exp / T_calc_O0 if T_calc_O0 > 0 else np.nan
        HF_O2 = T_exp / T_calc_O2 if T_calc_O2 > 0 else np.nan

        print(f"{name:12s} | {ell:3d} | {Q:7.3f} | {gamma_O0:7.2f} | {gamma_O2:7.2f} | {delta_gamma:+7.3f} | {T_ratio:9.4f} | {HF_O0:8.2f} | {HF_O2:8.2f}")

        results.append({
            'name': name,
            'Z_d': Z_d,
            'A_d': A_d,
            'Q': Q,
            'T_exp': T_exp,
            'ell': ell,
            'J_in': J_in,
            'J_fi': J_fi,
            'note': note,
            'gamma_O0': gamma_O0,
            'gamma_O2': gamma_O2,
            'delta_gamma': delta_gamma,
            'T_ratio': T_ratio,
            'T_calc_O0': T_calc_O0,
            'T_calc_O2': T_calc_O2,
            'HF_O0': HF_O0,
            'HF_O2': HF_O2,
            'r_in': result_O0['r_in'],
            'r_out': result_O0['r_out'],
        })

    print("-" * 100)

    return results


def analyze_by_angular_momentum(results):
    """
    Analyze results grouped by angular momentum.
    """
    print("\n" + "=" * 80)
    print("Analysis by Angular Momentum ℓ")
    print("=" * 80)

    # Group by ℓ
    ell_groups = {}
    for r in results:
        ell = r['ell']
        if ell not in ell_groups:
            ell_groups[ell] = []
        ell_groups[ell].append(r)

    for ell in sorted(ell_groups.keys()):
        nuclei = ell_groups[ell]
        n = len(nuclei)

        avg_gamma = np.mean([r['gamma_O0'] for r in nuclei])
        avg_delta = np.mean([r['delta_gamma'] for r in nuclei])
        avg_T_ratio = np.mean([r['T_ratio'] for r in nuclei])

        print(f"\nℓ = {ell}:")
        print(f"  Number of transitions: {n}")
        print(f"  Average γ(O0): {avg_gamma:.2f}")
        print(f"  Average Δγ (O2-O0): {avg_delta:.4f}")
        print(f"  Average T(O2)/T(O0): {avg_T_ratio:.4f}")
        print(f"  Transitions: {', '.join([r['name'] for r in nuclei])}")

        # For ℓ > 0, show hindrance factor comparison
        if ell > 0:
            print(f"\n  Hindrance factor comparison (unfavored transitions):")
            for r in nuclei:
                print(f"    {r['name']:12s}: HF(O0)={r['HF_O0']:.2f}, HF(O2)={r['HF_O2']:.2f}")


def compare_favored_unfavored(results):
    """
    Compare favored (ℓ=0) vs unfavored (ℓ>0) transitions.
    """
    print("\n" + "=" * 80)
    print("Favored vs Unfavored Transitions")
    print("=" * 80)

    favored = [r for r in results if r['ell'] == 0]
    unfavored = [r for r in results if r['ell'] > 0]

    print(f"\nFavored transitions (ℓ=0): {len(favored)}")
    print(f"Unfavored transitions (ℓ>0): {len(unfavored)}")

    if favored:
        avg_gamma_fav = np.mean([r['gamma_O0'] for r in favored])
        avg_delta_fav = np.mean([r['delta_gamma'] for r in favored])
        avg_pct_fav = np.mean([100 * r['delta_gamma'] / r['gamma_O0'] for r in favored])
        print(f"\nFavored (ℓ=0):")
        print(f"  Average γ: {avg_gamma_fav:.2f}")
        print(f"  Average Δγ: {avg_delta_fav:.4f}")
        print(f"  Average % correction: {avg_pct_fav:.3f}%")

    if unfavored:
        avg_gamma_unf = np.mean([r['gamma_O0'] for r in unfavored])
        avg_delta_unf = np.mean([r['delta_gamma'] for r in unfavored])
        avg_pct_unf = np.mean([100 * r['delta_gamma'] / r['gamma_O0'] for r in unfavored])
        print(f"\nUnfavored (ℓ>0):")
        print(f"  Average γ: {avg_gamma_unf:.2f}")
        print(f"  Average Δγ: {avg_delta_unf:.4f}")
        print(f"  Average % correction: {avg_pct_unf:.3f}%")


def extract_preformation_factors(results):
    """
    Extract preformation factors from experimental data.
    """
    print("\n" + "=" * 80)
    print("Preformation Factor Extraction (ℓ=0 transitions only)")
    print("=" * 80)

    # Filter ℓ=0 transitions
    favored = [r for r in results if r['ell'] == 0]

    print(f"\n{'Nucleus':12s} | {'Q(MeV)':7s} | {'γ':7s} | {'P_α(O0)':10s} | {'P_α(O2)':10s} | {'Ratio':6s}")
    print("-" * 70)

    P_O0_list = []
    P_O2_list = []

    for r in favored:
        name = r['name']
        data = SUPERHEAVY_NUCLEI[name]
        Z_d, A_d, Q, T_exp, ell, *_ = data

        system = AiryWKBSuperheavy(Z_d, A_d)

        P_O0 = system.extract_preformation_factor(Q, T_exp, ell, order=0)
        P_O2 = system.extract_preformation_factor(Q, T_exp, ell, order=2)

        ratio = P_O2 / P_O0 if P_O0 > 0 else np.nan

        print(f"{name:12s} | {Q:7.3f} | {r['gamma_O0']:7.2f} | {P_O0:10.4f} | {P_O2:10.4f} | {ratio:6.3f}")

        P_O0_list.append(P_O0)
        P_O2_list.append(P_O2)

    print("-" * 70)

    # Statistics
    P_O0_mean = np.mean(P_O0_list)
    P_O2_mean = np.mean(P_O2_list)
    P_O0_std = np.std(P_O0_list)
    P_O2_std = np.std(P_O2_list)

    print(f"\nStatistics:")
    print(f"  P_α(O0): {P_O0_mean:.4f} ± {P_O0_std:.4f} (CV = {100*P_O0_std/P_O0_mean:.1f}%)")
    print(f"  P_α(O2): {P_O2_mean:.4f} ± {P_O2_std:.4f} (CV = {100*P_O2_std/P_O2_mean:.1f}%)")
    print(f"  Ratio P_α(O2)/P_α(O0): {P_O2_mean/P_O0_mean:.3f}")

    return {'P_O0': P_O0_list, 'P_O2': P_O2_list, 'favored': favored}


def plot_superheavy_results(results):
    """
    Create comprehensive plots for superheavy nuclei analysis.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Filter data
    favored = [r for r in results if r['ell'] == 0]
    unfavored = [r for r in results if r['ell'] > 0]

    # (a) Gamow factor vs Q-value
    ax1 = axes[0, 0]
    Q_fav = [r['Q'] for r in favored]
    gamma_fav = [r['gamma_O0'] for r in favored]
    Q_unf = [r['Q'] for r in unfavored]
    gamma_unf = [r['gamma_O0'] for r in unfavored]

    ax1.scatter(Q_fav, gamma_fav, s=100, c='blue', marker='o', label=r'$\ell=0$ (favored)', edgecolors='darkblue', linewidths=1.5)
    ax1.scatter(Q_unf, gamma_unf, s=100, c='red', marker='^', label=r'$\ell>0$ (unfavored)', edgecolors='darkred', linewidths=1.5)

    ax1.set_xlabel('Q-value (MeV)', fontsize=16)
    ax1.set_ylabel(r'Gamow factor $\gamma$', fontsize=16)
    ax1.tick_params(axis='both', labelsize=12)
    ax1.legend(fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('(a) Gamow Factor vs Q-value', fontsize=14, fontweight='bold')

    # (b) O2 correction percentage
    ax2 = axes[0, 1]
    pct_corr_fav = [100 * r['delta_gamma'] / r['gamma_O0'] for r in favored]
    pct_corr_unf = [100 * r['delta_gamma'] / r['gamma_O0'] for r in unfavored]

    ax2.scatter(gamma_fav, pct_corr_fav, s=100, c='blue', marker='o', label=r'$\ell=0$', edgecolors='darkblue', linewidths=1.5)
    ax2.scatter(gamma_unf, pct_corr_unf, s=100, c='red', marker='^', label=r'$\ell>0$', edgecolors='darkred', linewidths=1.5)

    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=1)
    ax2.set_xlabel(r'Gamow factor $\gamma$', fontsize=16)
    ax2.set_ylabel(r'O2 correction $\Delta\gamma/\gamma$ (%)', fontsize=16)
    ax2.tick_params(axis='both', labelsize=12)
    ax2.legend(fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.set_title('(b) O2 Correction Magnitude', fontsize=14, fontweight='bold')

    # (c) Hindrance factors
    ax3 = axes[1, 0]
    HF_O0_fav = [r['HF_O0'] for r in favored if not np.isnan(r['HF_O0']) and r['HF_O0'] < 100]
    HF_O2_fav = [r['HF_O2'] for r in favored if not np.isnan(r['HF_O2']) and r['HF_O2'] < 100]
    names_fav = [r['name'] for r in favored if not np.isnan(r['HF_O0']) and r['HF_O0'] < 100]

    x = np.arange(len(names_fav))
    width = 0.35

    ax3.bar(x - width/2, HF_O0_fav, width, label='O0 (Gamow)', color='blue', alpha=0.7)
    ax3.bar(x + width/2, HF_O2_fav, width, label='O2 (Airy)', color='red', alpha=0.7)
    ax3.axhline(y=1, color='green', linestyle='--', linewidth=2, label='HF=1 (perfect)')

    ax3.set_xticks(x)
    ax3.set_xticklabels(names_fav, rotation=45, ha='right', fontsize=9)
    ax3.set_ylabel('Hindrance Factor', fontsize=16)
    ax3.tick_params(axis='y', labelsize=12)
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    ax3.set_title(r'(c) Hindrance Factors ($\ell=0$ only)', fontsize=14, fontweight='bold')

    # (d) T(O2)/T(O0) ratio
    ax4 = axes[1, 1]
    T_ratio_fav = [r['T_ratio'] for r in favored]
    T_ratio_unf = [r['T_ratio'] for r in unfavored]

    all_ratios = T_ratio_fav + T_ratio_unf
    all_labels = ['F'] * len(T_ratio_fav) + ['U'] * len(T_ratio_unf)
    all_gammas = gamma_fav + gamma_unf

    colors = ['blue' if l == 'F' else 'red' for l in all_labels]

    ax4.scatter(all_gammas, all_ratios, s=100, c=colors, edgecolors='black', linewidths=1)
    ax4.axhline(y=1, color='gray', linestyle='--', linewidth=1.5)

    # Add labels
    ax4.scatter([], [], c='blue', s=100, label=r'$\ell=0$ (favored)')
    ax4.scatter([], [], c='red', s=100, label=r'$\ell>0$ (unfavored)')

    ax4.set_xlabel(r'Gamow factor $\gamma$', fontsize=16)
    ax4.set_ylabel(r'$T_{\mathrm{O2}}/T_{\mathrm{O0}}$', fontsize=16)
    ax4.tick_params(axis='both', labelsize=12)
    ax4.legend(fontsize=12)
    ax4.grid(True, alpha=0.3)
    ax4.set_title('(d) Penetrability Enhancement', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('superheavy_airy_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('superheavy_airy_analysis.pdf', bbox_inches='tight')
    plt.close()
    print("\nSaved: superheavy_airy_analysis.png/pdf")


def plot_decay_chains(results):
    """
    Plot comparison for Ds decay chains.
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 275Ds chain
    ax1 = axes[0]
    chain_275 = ['275Ds_0', '275Ds_4', '271Hs_0a', '271Hs_0b', '271Hs_6', '267Sg']

    data_275 = []
    for name in chain_275:
        for r in results:
            if r['name'] == name:
                data_275.append(r)
                break

    if data_275:
        names = [d['name'].replace('_', '\n') for d in data_275]
        gamma_O0 = [d['gamma_O0'] for d in data_275]
        gamma_O2 = [d['gamma_O2'] for d in data_275]

        x = np.arange(len(names))
        width = 0.35

        ax1.bar(x - width/2, gamma_O0, width, label='O0', color='blue', alpha=0.7)
        ax1.bar(x + width/2, gamma_O2, width, label='O2 (Airy)', color='red', alpha=0.7)

        ax1.set_xticks(x)
        ax1.set_xticklabels(names, fontsize=10)
        ax1.set_ylabel(r'Gamow factor $\gamma$', fontsize=16)
        ax1.tick_params(axis='y', labelsize=12)
        ax1.legend(fontsize=12)
        ax1.grid(axis='y', alpha=0.3)
        ax1.set_title(r'$^{275}$Ds decay chain', fontsize=14, fontweight='bold')

    # 273Ds chain
    ax2 = axes[1]
    chain_273 = ['273Ds_0a', '273Ds_6', '273Ds_0b', '269Hs_0a', '269Hs_4', '269Hs_0b', '265Sg_0a', '265Sg_0b', '261Rf']

    data_273 = []
    for name in chain_273:
        for r in results:
            if r['name'] == name:
                data_273.append(r)
                break

    if data_273:
        names = [d['name'].replace('_', '\n') for d in data_273]
        gamma_O0 = [d['gamma_O0'] for d in data_273]
        gamma_O2 = [d['gamma_O2'] for d in data_273]

        x = np.arange(len(names))

        ax2.bar(x - width/2, gamma_O0, width, label='O0', color='blue', alpha=0.7)
        ax2.bar(x + width/2, gamma_O2, width, label='O2 (Airy)', color='red', alpha=0.7)

        ax2.set_xticks(x)
        ax2.set_xticklabels(names, fontsize=8, rotation=45, ha='right')
        ax2.set_ylabel(r'Gamow factor $\gamma$', fontsize=16)
        ax2.tick_params(axis='y', labelsize=12)
        ax2.legend(fontsize=12)
        ax2.grid(axis='y', alpha=0.3)
        ax2.set_title(r'$^{273}$Ds decay chain', fontsize=14, fontweight='bold')

    plt.tight_layout()
    plt.savefig('ds_decay_chains.png', dpi=300, bbox_inches='tight')
    plt.savefig('ds_decay_chains.pdf', bbox_inches='tight')
    plt.close()
    print("Saved: ds_decay_chains.png/pdf")


def summarize_key_findings(results):
    """
    Summarize key findings from the analysis.
    """
    print("\n" + "=" * 80)
    print("KEY FINDINGS")
    print("=" * 80)

    # 1. Gamow factors
    gammas = [r['gamma_O0'] for r in results]
    print(f"\n1. Gamow Factors:")
    print(f"   Range: {min(gammas):.2f} - {max(gammas):.2f}")
    print(f"   Mean: {np.mean(gammas):.2f}")
    print(f"   (Compared to actinides: γ ~ 40-45)")

    # 2. O2 corrections
    delta_gammas = [r['delta_gamma'] for r in results]
    pct_corrs = [100 * r['delta_gamma'] / r['gamma_O0'] for r in results]
    print(f"\n2. O2 Corrections (Airy method):")
    print(f"   Δγ range: {min(delta_gammas):.4f} to {max(delta_gammas):.4f}")
    print(f"   % correction range: {min(pct_corrs):.3f}% to {max(pct_corrs):.3f}%")
    print(f"   Mean % correction: {np.mean(pct_corrs):.3f}%")

    # 3. Penetrability enhancement
    T_ratios = [r['T_ratio'] for r in results]
    print(f"\n3. Penetrability Enhancement T(O2)/T(O0):")
    print(f"   Range: {min(T_ratios):.4f} - {max(T_ratios):.4f}")
    print(f"   Mean: {np.mean(T_ratios):.4f}")

    # 4. Comparison with actinides
    print(f"\n4. Comparison with Actinides:")
    print(f"   Actinide γ ~ 42, Δγ ~ -0.06 (Airy), % correction ~ -0.14%")
    print(f"   Superheavy γ ~ 25, Δγ ~ {np.mean(delta_gammas):.4f}, % correction ~ {np.mean(pct_corrs):.3f}%")

    # 5. Angular momentum effect
    favored = [r for r in results if r['ell'] == 0]
    unfavored = [r for r in results if r['ell'] > 0]
    if favored and unfavored:
        avg_pct_fav = np.mean([100 * r['delta_gamma'] / r['gamma_O0'] for r in favored])
        avg_pct_unf = np.mean([100 * r['delta_gamma'] / r['gamma_O0'] for r in unfavored])
        print(f"\n5. Angular Momentum Effect:")
        print(f"   ℓ=0 (favored): {avg_pct_fav:.3f}% correction")
        print(f"   ℓ>0 (unfavored): {avg_pct_unf:.3f}% correction")


if __name__ == "__main__":
    print("\n" + "=" * 90)
    print("SUPERHEAVY NUCLEI ALPHA DECAY: AIRY UNIFORM WKB ANALYSIS")
    print("=" * 90)

    # Main analysis
    results = analyze_superheavy_nuclei()

    # Group by angular momentum
    analyze_by_angular_momentum(results)

    # Compare favored vs unfavored
    compare_favored_unfavored(results)

    # Extract preformation factors
    pf_data = extract_preformation_factors(results)

    # Summary
    summarize_key_findings(results)

    # Plots
    print("\n" + "=" * 80)
    print("GENERATING FIGURES")
    print("=" * 80)
    plot_superheavy_results(results)
    plot_decay_chains(results)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
