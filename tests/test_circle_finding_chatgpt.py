import numpy as np


# =====================================================
# FORCE METHOD
# =====================================================

def force_circle_fit(points, steps):
    center = np.mean(points, axis=0)
    radius = np.mean(np.linalg.norm(points - center, axis=1))
    last_forces = np.zeros(2)
    n = len(points)

    for _ in range(steps):
        distances = np.linalg.norm(points - center, axis=1)
        forces = np.zeros(2)

        for p, d in zip(points, distances):
            forces += (p - center) * (d - radius) / (d + 1e-9)

        forces /= n

        dot = np.dot(forces, last_forces)
        #denom = (np.linalg.norm(forces) * np.linalg.norm(last_forces) + 1e-9)
        #damp_factor = (dot / denom) > 0
        damp_factor = dot > 0
        last_forces *= damp_factor
        last_forces += forces

        center += last_forces
        radius = np.mean(np.linalg.norm(points - center, axis=1))

    return center, radius


# =====================================================
# KÅSA + GAUSS–NEWTON
# =====================================================

def fit_circle_kasa(points):
    x = points[:,0]
    y = points[:,1]
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x*x + y*y)
    params, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    a, b, c = params
    center = np.array([-a/2, -b/2])
    radius = np.sqrt((a*a + b*b)/4 - c)
    return center, radius

def refine_gauss_newton(points, center, radius, iters):
    cx, cy = center
    r = radius
    for _ in range(iters):
        J = []
        residuals = []
        for px, py in points:
            dx = px - cx
            dy = py - cy
            dist = np.sqrt(dx*dx + dy*dy) + 1e-9
            res = dist - r
            J.append([-dx/dist, -dy/dist, -1])
            residuals.append(res)
        J = np.array(J)
        residuals = np.array(residuals)
        delta, _, _, _ = np.linalg.lstsq(J, -residuals, rcond=None)
        cx += delta[0]
        cy += delta[1]
        r  += delta[2]
    return np.array([cx, cy]), r


# =====================================================
# METRICS
# =====================================================

def metrics(points, center, radius, true_center, true_radius):
    center_err = np.linalg.norm(center - true_center)
    radius_err = abs(radius - true_radius)
    residuals = np.linalg.norm(points - center, axis=1) - radius
    rmse = np.sqrt(np.mean(residuals**2))
    return center_err, radius_err, rmse


# =====================================================
# SAMPLING: NOISY + PARTIAL ARCS
# =====================================================

def sample_circle(n=30, noise_sigma=0.25, arc_fraction=0.25):
    """
    arc_fraction: 1.0 = full circle
                  0.25 = 90° arc
                  0.10 = 36° arc
    """
    true_center = np.random.randn(2) * 2
    true_radius = np.random.rand()*2 + 2

    arc = 2*np.pi * arc_fraction
    start = np.random.uniform(0, 2*np.pi)
    angles = np.random.uniform(start, start+arc, n)

    points = true_center + np.column_stack([np.cos(angles), np.sin(angles)]) * true_radius
    points += np.random.randn(n, 2) * noise_sigma

    return points, true_center, true_radius


# =====================================================
# EXPERIMENT TABLE
# =====================================================

def run_table(
    trials=300,
    force_steps_list=[1, 2, 3, 5, 10, 20],
    gn_steps_list=[0, 1, 2],
    n_points=20,
    noise_sigma=0.25,
    arc_fraction=0.25,
):

    print("\n===============================================================")
    print("                NOISY / PARTIAL-ARC CIRCLE FIT")
    print(f"         {noise_sigma=}, {arc_fraction=}  ({int(360*arc_fraction)}° arc)")
    print(f"         averaged over {trials} trials, n_points={n_points}")
    print("===============================================================\n")

    print(f"{'Method':<20} {'Steps':<6}  {'CenterErr (mean±std)':<26} {'RadiusErr (mean±std)':<26} {'RMSE (mean±std)'}")
    print("-"*100)

    # -------------------------------------------------
    # Force method
    # -------------------------------------------------
    for fs in force_steps_list:
        ce, re, rm = [], [], []
        for _ in range(trials):
            points, tc, tr = sample_circle(n=n_points, noise_sigma=noise_sigma, arc_fraction=arc_fraction)
            c, r = force_circle_fit(points, fs)
            a, b, c_rm = metrics(points, c, r, tc, tr)
            ce.append(a); re.append(b); rm.append(c_rm)

        print(f"{'Force':<20} {fs:<6}  "
              f"{np.mean(ce):.4f}±{np.std(ce):.4f}      "
              f"{np.mean(re):.4f}±{np.std(re):.4f}      "
              f"{np.mean(rm):.4f}±{np.std(rm):.4f}")

    # -------------------------------------------------
    # Kasa and GN
    # -------------------------------------------------
    for gs in gn_steps_list:
        ce, re, rm = [], [], []
        for _ in range(trials):
            points, tc, tr = sample_circle(n=n_points, noise_sigma=noise_sigma, arc_fraction=arc_fraction)
            center, radius = fit_circle_kasa(points)
            if gs > 0:
                center, radius = refine_gauss_newton(points, center, radius, gs)
            a, b, c_rm = metrics(points, center, radius, tc, tr)
            ce.append(a); re.append(b); rm.append(c_rm)

        method_name = "Kasa" if gs == 0 else "Kasa+GN"
        print(f"{method_name:<20} {gs:<6}  "
              f"{np.mean(ce):.4f}±{np.std(ce):.4f}      "
              f"{np.mean(re):.4f}±{np.std(re):.4f}      "
              f"{np.mean(rm):.4f}±{np.std(rm):.4f}")


if __name__ == "__main__":
    run_table()
