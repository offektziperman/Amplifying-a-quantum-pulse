import numpy as np
from idesolver import IDESolver

OMEGA_MAX = 100
N_omega = 1000
Z_MAX = 10
N_z = 100


def find_UV(init_guess, f):
    pass

def get_pump(width, amplitude, omegas):
    OMEGA_1, OMEGA_2 = np.meshgrid(omegas, omegas)
    OMEGA_SUM = OMEGA_1 + OMEGA_2
    pump = amplitude / np.sqrt(2 * np.pi * width ** 2) * np.exp(-OMEGA_SUM ** 2 / (2 * width ** 2))

def single_cycle(U, V, f):
    solver = IDESolver(
        omegas=np.linspace(-OMEGA_MAX, OMEGA_MAX, N_omega),
        y_0=0,
        c=lambda x, y: y - (.5 * x) + (1 / (1 + x)) - np.log(1 + x),
        d=lambda x: 1 / (np.log(2)) ** 2,
        k=lambda x, s: x / (1 + s),
        f=lambda y: y,
        lower_bound=lambda x: 0,
        upper_bound=lambda x: 1,
    )


import numpy as np
import matplotlib.pyplot as plt
from idesolver import IDESolver

c = 1
z = np.linspace(0, Z_MAX, N_z)
omegas_1 = np.linspace(-OMEGA_MAX, OMEGA_MAX, N_omega)
omegas_2 = np.linspace(-OMEGA_MAX, OMEGA_MAX, N_omega)
k0 = np.zeros([len(omegas_1), len(omegas_2)])

U_0 = np.zeros([len(omegas_1), len(omegas_2)])
V_0 = np.zeros([len(omegas_1), len(omegas_2)])
U_0_vec, V_0_vec = U_0.flatten(), V_0.flatten()
OMEGAS_1, OMEGAS_2 = np.meshgrid(omegas_1, omegas_2)
delta_k = OMEGAS_1 / c - OMEGAS_2 / c



pump = get_pump(1, 1, omegas_1)


def k(z, omega_p):



y0 = np.concatenate([U_0_vec, V_0_vec])
solver = IDESolver(
    x=z,
    y_0=y0,
    c=lambda x, y: 0,
    d=lambda x: 1,
    k=k,
    f=lambda y: y,
    lower_bound=lambda x: 0,
    upper_bound=lambda x: x
)

solver.solve()

fig = plt.figure(dpi=100)
ax = fig.add_subplot(111)

exact = [np.sin(solver.x), np.cos(solver.x)]

ax.plot(solver.x, solver.y[1], label='IDESolver Solution', linestyle='-', linewidth=3)
ax.plot(solver.x, exact[1], label='Analytic Solution', linestyle=':', linewidth=3)

ax.legend(loc='best')
ax.grid(True)

ax.set_title(f'Solution for Global Error Tolerance = {solver.global_error_tolerance}')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y(x)$')

plt.show()
