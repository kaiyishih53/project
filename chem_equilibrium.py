from dataclasses import dataclass
import math


@dataclass
class WeakAcidSystem:
    """
    弱酸解離系統：
    HA ⇌ H+ + A-

    參數：
    - c0: 初始 HA 濃度 (mol/L)
    - Ka: 酸解離常數
    """
    c0: float   # 初始 HA 濃度
    Ka: float   # 酸解離常數

    def solve_equilibrium(self):
        """
        利用精確二次方程式解：
        Ka = x^2 / (c0 - x)
        => Ka(c0 - x) = x^2
        => Ka*c0 - Ka*x - x^2 = 0
        => x^2 + Ka*x - Ka*c0 = 0
        """
        a = 1.0
        b = self.Ka
        c = - self.Ka * self.c0

        disc = b * b - 4 * a * c
        if disc < 0:
            raise ValueError("無實數解，請檢查輸入參數。")

        x1 = (-b + math.sqrt(disc)) / (2 * a)
        x2 = (-b - math.sqrt(disc)) / (2 * a)

        # 物理解要求：0 ≤ x ≤ c0
        candidates = [x for x in (x1, x2) if 0 <= x <= self.c0]
        if not candidates:
            raise ValueError("沒有物理上合理的解，請檢查濃度與 Ka。")

        x = max(candidates)  # 通常較大的那個根是合理解

        h_conc = x          # [H+]
        a_minus = x         # [A-]
        ha = self.c0 - x    # [HA]

        return {
            "[H+]": h_conc,
            "[A-]": a_minus,
            "[HA]": ha,
            "pH": -math.log10(h_conc)
        }


def run_weak_acid_mode():
    print("=== 弱酸解離平衡計算器 (HA ⇌ H+ + A-) ===")
    try:
        c0 = float(input("請輸入初始酸濃度 c0 (mol/L): "))
        Ka = float(input("請輸入酸解離常數 Ka: "))
    except ValueError:
        print("輸入格式錯誤，請輸入數字。")
        return

    if c0 <= 0 or Ka <= 0:
        print("c0 與 Ka 必須為正數。")
        return

    system = WeakAcidSystem(c0=c0, Ka=Ka)
    try:
        result = system.solve_equilibrium()
    except ValueError as e:
        print("計算失敗：", e)
        return

    print("\n--- 計算結果 ---")
    print(f"[HA]  = {result['[HA]']:.6e} mol/L")
    print(f"[A-]  = {result['[A-]']:.6e} mol/L")
    print(f"[H+]  = {result['[H+]']:.6e} mol/L")
    print(f"pH    = {result['pH']:.3f}")
    print("----------------\n")


def main():
    while True:
        print("===== 化學平衡計算器 =====")
        print("1) 弱酸解離平衡 (HA ⇌ H+ + A-)")
        # 未來可以在這裡加：
        # 2) 多元酸平衡
        # 3) 鹽類水解
        # 4) 溶解度積 Ksp 平衡
        print("0) 離開")
        choice = input("請選擇功能編號: ").strip()

        if choice == "1":
            run_weak_acid_mode()
        elif choice == "0":
            print("感謝使用，掰掰～")
            break
        else:
            print("無效選項，請重新輸入。\n")


if __name__ == "__main__":
    main()
