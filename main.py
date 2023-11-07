import matplotlib.pyplot as plt
import numpy as np
import qutip
from scipy.integrate import quad
from cmath import phase
from numpy import arccosh, abs, arccos, conj, sin
from qutip import destroy, qeye, tensor, squeeze, basis, ket2dm, displace
from scipy.interpolate import interp1d
import scipy
import pickle
import matplotlib.colors as mcolors

from calc_PDC import calc_PDC

tp = 4
tau = 1
gamma = 1
t0 = 10
T = 30
N_se = 50
N = 50

e = lambda t: np.exp(t)
sinh = lambda t: np.sinh(t)
cosh = lambda t: np.cosh(t)
g = lambda t: 1


def find_closest_index_sorted(arr, target):
    left = 0
    right = len(arr) - 1

    while left <= right:
        mid = (left + right) // 2

        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    # At this point, left > right
    if right < 0:
        return left
    elif left >= len(arr):
        return right
    else:
        if abs(arr[right] - target) < abs(arr[left] - target):
            return right
        else:
            return left


def g1_t(times, epsilon, u, rho_u):
    '''

    :return:
    '''
    N = int(np.size(rho_u) ** 0.5)
    a = destroy(N)
    corr_n = qutip.expect(a.dag() * a, rho_u)
    corr_aa = qutip.expect(a * a, rho_u)
    corr_a_dag_a_dag = qutip.expect(a.dag() * a.dag(), rho_u)

    dt = times[1] - times[0]
    Nt = len(times)
    g1 = np.zeros([Nt, Nt], dtype=complex)
    g1_input_dependant = np.zeros([Nt, Nt], dtype=complex)

    C = lambda t1, t2: np.sqrt(gamma) * cosh(quad(epsilon, t1, t2
                                                  )[0]) * e(1 / 2 * gamma * (t2 - t1))
    D = lambda t1, t2: np.sqrt(gamma) * sinh(quad(epsilon, t1, t2)[0]) * e(1 / 2 * gamma * (t2 - t1))
    C_array = np.array([[C(t1, t2) for t1 in times] for t2 in times]).T
    D_array = np.array([[D(t1, t2) for t1 in times] for t2 in times]).T
    D_func_array = np.array([
        [np.sum([conj(D_array[i, k]) * D_array[j, k] * dt for k in range(min(i, j))]) for i in range(len(times))] for j
        in
        range(len(times))])
    C_u_array = np.array([np.sum([conj(C_array[i, k] * u[k]) * dt for k in range(i)]) for i in range(len(times))])
    D_u_array = np.array([np.sum([conj(D_array[i, k] * u[k]) * dt for k in range(i)]) for i in range(len(times))])

    for i, t1 in enumerate(times):
        for j, t2 in enumerate(times):
            g1[i, j] = conj(D_array[i, 0]) * D_array[j, 0] + D_func_array[i, j] + \
                       corr_n * (conj(u[i]) * u[j] - conj(u[i]) * C_u_array[j] - u[j] * conj(C_u_array[i]) \
                                 + conj(C_u_array[i]) * C_u_array[j] + conj(D_u_array[i]) * D_u_array[j]) \
                       - corr_aa * conj(u[i]) * D_u_array[j] - corr_a_dag_a_dag * u[j] * conj(D_u_array[i]) \
                       + corr_aa * C_u_array[i] * conj(D_u_array[j]) + \
                       corr_a_dag_a_dag * D_u_array[i] * conj(C_u_array[j])
            g1_input_dependant[i, j] = corr_n * (conj(u[i]) * u[j] + (
                    -conj(u[i]) * C_u_array[j] - u[j] * conj(C_u_array[i])) + \
                                                 conj(C_u_array[i]) * C_u_array[j] \
                                                 + conj(D_u_array[i]) * D_u_array[j]) - \
                                       corr_aa * conj(u[i]) * D_u_array[j] - corr_a_dag_a_dag * u[j] * conj(
                D_u_array[i]) \
                                       + corr_aa * C_u_array[i] * conj(D_u_array[j]) + \
                                       corr_a_dag_a_dag * D_u_array[i] * conj(C_u_array[j])
    return g1, g1_input_dependant


def get_occupation_of_modes(g1_mat, num_modes=2):
    '''

    :return:
    '''
    w, v = np.linalg.eig(g1_mat)  # diagonalize gamma matrix
    occupations = np.sort(np.abs(w))[-1:-num_modes - 1:-1]
    return occupations


def get_most_occupied_mode_function(times, g1_mat, num_modes=1):
    '''

    :return:
    '''
    dt = times[1] - times[0]
    w, v = np.linalg.eig(g1_mat)  # diagonalize gamma matrix
    # plt.scatter(np.linspace(1, len(occupations), len(occupations)), occupations)
    # plt.show()
    mode_1 = v[:, np.abs(w) == sorted(np.abs(w), reverse=True)[0]]
    mode_1 = np.real(mode_1) / np.sqrt(sum(np.abs(mode_1) ** 2 * dt))

    mode_2 = v[:, np.abs(w) == sorted(np.abs(w), reverse=True)[1]]
    mode_2 = np.real(mode_2) / np.sqrt(sum(np.abs(mode_2) ** 2 * dt))
    mode = interp1d(times, mode_1.T)
    mode_func1 = lambda t: mode(t)[0] if t < np.max(times) and t > 0 else 0

    if num_modes == 1:
        return mode_func1
    if num_modes == 2:
        return mode_1, mode_2


def find_coefficiants_for_transformation_arrays(times, u, v_func, epsilon):
    '''

    :param u:
    :param v:
    :param epsilon_int:
    :return:
    '''
    v = np.array([v_func(t) for t in times])
    dt = times[2] - times[1]

    epsilon_int_ar = np.cumsum([epsilon(t) * dt for t in times])

    def epsilon_int(t1, t2):
        index_t1 = find_closest_index_sorted(times, t1)
        index_t2 = find_closest_index_sorted(times, t2)
        return epsilon_int_ar[index_t2] - epsilon_int_ar[index_t1]

    f = v.conjugate() * gamma * np.exp(-gamma * times / 2) * cosh(epsilon_int_ar)
    dt = times[1] - times[0]

    A = np.sum(f) * dt
    f = v.conjugate() * gamma * np.exp(-gamma * times / 2) * sinh(epsilon_int_ar)
    B = np.sum(f) * dt
    second_term = np.zeros_like(times)

    for i, tprime in enumerate(times):
        second_term[i] = e(gamma / 2 * tprime) * sum(np.array([v[j].conjugate() * gamma * e(-gamma / 2 * t)
                                                               * cosh(epsilon_int(tprime, t)) for j, t in
                                                               enumerate(times)])[times > tprime]) * dt

    plt.show()
    c_tilde = v.conjugate() - second_term

    # plt.plot(times, u)
    # plt.plot(times, c_tilde)
    # plt.show()
    C = (sum(np.abs(c_tilde) ** 2) * dt) ** 0.5

    d_tilde = np.zeros_like(times)
    for i, tprime in enumerate(times):
        d_tilde[i] = - sum(np.array([v[j].conjugate() * gamma * np.exp(gamma / 2 * (tprime - t))
                                     * sinh(epsilon_int(tprime, t)) for j, t in enumerate(times)])[times > tprime]) * dt
    D = (sum(np.abs(d_tilde) ** 2) * dt) ** 0.5

    c = c_tilde / C
    d = d_tilde / D

    # overlap = lambda f, g: quad(lambda t: f(t).conjugate() * g(t), 0, np.inf)[0]

    overlap_arr = lambda f, g: sum(conj(f) * g * dt)
    cu, ud = overlap_arr(c, u), overlap_arr(u, d)

    h = (c - u * cu.conjugate()) / np.sqrt(1 - cu * cu.conjugate())
    k = (d - u * ud) / np.sqrt(1 - ud * ud.conjugate())

    kh = overlap_arr(k, h)

    coef_array = np.asanyarray([[C * cu, D * ud], [C * np.sqrt(1 - abs(cu) ** 2) * kh, D * np.sqrt(1 - abs(ud) ** 2)],
                                [C * np.sqrt(1 - abs(cu) ** 2) * np.sqrt(1 - abs(kh) ** 2)]], dtype=object)
    return coef_array


def squeezing_and_bs_from_coefs(A, B):
    '''

    :return:
    '''
    r = arccosh(abs(A) / np.sqrt(abs(A) ** 2 - abs(B) ** 2))
    return r, phase(B) - phase(A), phase(A), np.sqrt(abs(A) ** 2 - abs(B) ** 2)


def rotate_density_mat(rho, phi):
    '''

    :param rho:
    :param phi:
    :return:
    '''
    a = destroy(N)
    n = a.dag() * a
    rot = (-1j * n * phi).expm()
    return rot * rho * rot.dag()


def mix_modes(rho_a, rho_b, a_coef, b_coef):
    '''

    :param rho_a:
    :param rho_b:
    :param a_coef:
    :param b_coef:
    :return:
    '''
    N = np.size(rho_b, 1)
    a = tensor(destroy(N), qeye(N))
    b = tensor(qeye(N), destroy(N))
    rho_a, rho_b = rotate_density_mat(rho_a, phase(a_coef)), rotate_density_mat(rho_b, phase(b_coef))
    norm_a_coef = np.abs(a_coef) / np.sqrt(abs(a_coef) ** 2 + abs(b_coef) ** 2)
    theta = arccos(norm_a_coef)
    to_exp = -1j * theta * (a.dag() * b + a * b.dag())
    U_bs = to_exp.expm()
    rho_t = tensor(rho_a, rho_b)
    rho_t = U_bs * rho_t * U_bs.dag()
    return rho_t.ptrace(0), np.sqrt(abs(a_coef) ** 2 + abs(b_coef) ** 2)


def bloch_messiah_transform(rho_A, rho_B, theta_1, phi_1, r_1, r_2, theta_2, phi_2):
    rho_t = qutip.tensor(rho_A, rho_B)
    N = np.size(rho_B, 1)
    a = tensor(destroy(N), qeye(N))
    b = tensor(qeye(N), destroy(N))
    to_exp = -1j * theta_1 * (a.dag() * b * np.exp(1j * phi_1) + a * b.dag() * np.exp(-1j * phi_1))
    U_bs = to_exp.expm()
    rho_t = U_bs * rho_t * U_bs.dag()
    S1 = tensor(squeeze(N, r_1), qeye(N))
    S2 = tensor(qeye(N), squeeze(N, r_2))

    rho_t = S2 * S1 * rho_t * S1.dag() * S2.dag()
    to_exp = -1j * theta_2 * (a.dag() * b * np.exp(1j * phi_2) + a * b.dag() * np.exp(-1j * phi_2))
    U_bs = to_exp.expm()
    rho_t = U_bs * rho_t * U_bs.dag()
    return rho_t.ptrace(0)


def is_unitary(matrix):
    conjugate_transpose = np.conj(matrix.T)
    product = np.dot(matrix, conjugate_transpose)
    identity = np.eye(matrix.shape[0])

    return np.allclose(product, identity)


def find_bloch_messiah_reduction(A, B, C, D):
    w1 = 1 / (abs(A) ** 2 + abs(C) ** 2) ** 0.5 * np.array([[A, -C], [C, A]], dtype=complex)
    w2 = 1 / (abs(B) ** 2 + abs(D) ** 2) ** 0.5 * np.array([[B, D], [D, -B]], dtype=complex)

    mat_for_E = conj(w1.T) @ conj(w2)
    E = scipy.linalg.sqrtm(mat_for_E)

    sigma_a = np.array([[abs(A) ** 2 + abs(C) ** 2, 0], [0, abs(A) ** 2 + abs(C) ** 2]], dtype=complex) ** 0.5
    sigma_b = np.array([[abs(B) ** 2 + abs(D) ** 2, 0], [0, abs(B) ** 2 + abs(D) ** 2]], dtype=complex) ** 0.5

    return E, [sigma_a, sigma_b], conj(w2 @ E)


def straighten_phase_cat(rho, light_dim):
    N = np.size(rho[0, :])
    psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
    psi_i = psi_i / psi_i.norm()
    rho_cat = qutip.ket2dm(psi_i)

    fidelity_min = qutip.fidelity(qutip.Qobj(rho_cat), qutip.Qobj(rho))
    a = qutip.destroy(light_dim)
    n = a.dag() * a
    phi_0_max = 0
    rho_max = rho
    for phi_0 in np.linspace(0, 2 * np.pi, 500):
        rho_rot = (1j * phi_0 * n).expm() * rho * (-1j * phi_0 * n).expm()
        fidelity_cat = qutip.fidelity(qutip.Qobj(rho_cat), qutip.Qobj(rho))
        if fidelity_cat > fidelity_min:
            fidelity_min = fidelity_cat
            rho_max = rho_rot
    return rho_max


def straighten_phase(rho: qutip.Qobj, light_dim: int, squeeze: bool = False) -> qutip.Qobj:
    '''

    :param squeeze:
    :param rho:
    :param light_dim:
    :return:
    '''
    psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
    psi_i = psi_i / psi_i.norm()
    rho_cat = qutip.ket2dm(psi_i)
    a = qutip.destroy(light_dim)
    x = a + a.dag()
    x2_expect_min = np.trace(x * x * rho)
    a = qutip.destroy(light_dim)
    n = a.dag() * a
    phi_0_max = 0
    rho_max = rho
    for phi_0 in np.linspace(0, 2 * np.pi, 100):
        rho_rot = (1j * phi_0 * n).expm() * rho * (-1j * phi_0 * n).expm()
        x2_expect = np.trace(x * x * rho_rot) - np.trace(x * rho_rot) ** 2
        if x2_expect < x2_expect_min:
            x2_expect_min = x2_expect
            rho_max = rho_rot
    return rho_max


def get_params_from_reduction(U, sigma_a, V):
    '''

    :param U:
    :param sigma_a:
    :param V:
    :return:
    '''

    theta_2 = arccos(np.abs(U[0, 0]))
    phi_2 = phase(1j * U[0, 1]) - phase(U[0, 0])

    r_1 = arccosh(sigma_a[0, 0])
    r_2 = arccosh(sigma_a[1, 1])
    V_dag = conj(V.T)
    theta_1 = arccos(np.abs(V[0, 0]))
    phi_1 = phase(1j * V_dag[0, 1]) - phase(V_dag[1, 1])

    return theta_1, phi_1, r_1, r_2, theta_2, phi_2


def transform_quantum_state(rho_array, coef_array, THRESH=0.02):
    '''

    :param rho_array:
    :param coef_array:
    :return:
    '''
    A, B, C, D = coef_array[0][0], coef_array[0][1], coef_array[1][0], coef_array[1][1]
    z = (np.abs(A) ** 2 + np.abs(C) ** 2 - np.abs(B) ** 2 - np.abs(D) ** 2) ** 0.5
    A, B, C, D = (A, B, C, D) / z

    U, [sigma_a, sigma_b], V = find_bloch_messiah_reduction(A, B, C, D)

    theta_1, phi_1, r_1, r_2, theta_2, phi_2 = get_params_from_reduction(U, sigma_a, V)
    rho = bloch_messiah_transform(rho_array[0], rho_array[1], theta_1, phi_1, r_1, r_2, theta_2, phi_2)
    rho, coef = mix_modes(rho, rho_array[2], z, coef_array[2][0])
    return rho


def plot_figure_4():
    f, axes = plt.subplots(3, 5, sharex=True, sharey=True, figsize=(4, 2.6))
    sigma_e = np.linspace(0.3, 10, 50)

    purity = np.zeros([3, len(sigma_e)])
    population = np.zeros([3, len(sigma_e)])
    plt.subplots_adjust(wspace=0.03, hspace=0.03)

    s = [0, 4, 10, 49]
    s = [0,4]
    sigma_e = sigma_e[s]
    print(sigma_e)
    for i, st in enumerate(['Vacuum', 'Fock', 'Cat']):
        for j, sig in enumerate(sigma_e):
            rho = np.load('density_mats/' + str(st) + '/eps_index=' + str(s[j]) + '.npy')
            # if i == 2 and j == 4:
            N_r = np.size(rho[0, :])
            if i == 2:
                rho = straighten_phase_cat(rho, N_r)
            plt_wgnr(rho, f, axes[i, j])

            axes[i, j].set_xlim([-3.5, 3.5])
            axes[i, j].set_ylim([-6, 6])

            axes[i, j].set_aspect('equal')
    # for i, st in enumerate(['Vacuum', 'Fock', 'Cat']):
    #     rho = np.load('density_mats/' + str(st) + '/times_5.npy')
    #     N_r = np.size(rho[0, :])
    #     plt_wgnr(rho, f, axes[i, 4])
    #
    #     axes[i, 4].set_xlim([-3.5, 3.5])
    #     axes[i, 4].set_ylim([-6, 6])
    #     axes[i, 4].set_aspect('equal')

    f.set_figwidth(5)
    f.set_figheight(5)
    plt.savefig('Wigners.svg')
    plt.show()


def find_fidelity_and_squeezing(rho, rho_i):
    N = np.size(rho[0, :])
    max_xi = 0
    max_fid_to_squeezed = 0
    for xi in np.linspace(0, 2, 400):
        rho_squeezed = squeeze(N, xi) * qutip.Qobj(rho_i) * squeeze(N, xi).dag()

        fid_to_squeezed = qutip.fidelity(rho_squeezed, qutip.Qobj(rho))
        if fid_to_squeezed > max_fid_to_squeezed:
            max_fid_to_squeezed = fid_to_squeezed
            max_xi = xi
    return max_fid_to_squeezed, max_xi


def times_5_gain():
    N_se = 50
    rho_u = qutip.ket2dm(qutip.basis(N, 0))
    w = 10
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))

    A = quad(epsilon, 0, T)[0]
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5 * 5
    rho = calculate_final_density_matrix(epsilon, rho_u)
    np.save('density_mats/Vacuum/times_5', rho)
    rho_u = qutip.ket2dm(qutip.basis(N, 1))
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))
    A = quad(epsilon, 0, T)[0]
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5 * 5
    rho = calculate_final_density_matrix(epsilon, rho_u)
    np.save('density_mats/Fock/times_5', rho)

    psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
    psi_i = psi_i / psi_i.norm()
    rho_u = qutip.ket2dm(psi_i)
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))
    A = quad(epsilon, 0, T)[0]
    epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5 * 5
    rho = calculate_final_density_matrix(epsilon, rho_u)
    np.save('density_mats/Cat/times_5', rho)


def calculate_figure_4_lineplots():
    sigma_e = np.linspace(0.3, 10, 50)

    purity = np.zeros([3, len(sigma_e)])
    population = np.zeros([3, len(sigma_e)])
    squeezing = np.zeros([3, len(sigma_e)])
    fidelity_to_squeezed = np.zeros([3, len(sigma_e)])
    for i, st in enumerate(['Vacuum', 'Fock', 'Cat']):
        if i == 0:
            rho = np.load('density_mats/Vacuum/eps_index=0.npy')
            N = np.size(rho[0, :])
            rho_i = qutip.ket2dm(basis(N, 0))
        if i == 1:
            rho = np.load('density_mats/Fock/eps_index=0.npy')
            N = np.size(rho[0, :])
            rho_i = qutip.ket2dm(basis(N, 1))
        if i == 2:
            rho = np.load('density_mats/Cat/eps_index=0.npy')
            N = np.size(rho[0, :])
            psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
            psi_i = psi_i / psi_i.norm()
            rho_i = qutip.ket2dm(psi_i)

        for j, sig in enumerate(sigma_e):
            if i == 2 and j == 26:
                continue

            print(j)
            rho = np.load('density_mats/' + str(st) + '/eps_index=' + str(j) + '.npy')
            N = np.size(rho, 1)
            a = qutip.destroy(N)
            n = a.dag() * a

            purity[i, j] = np.trace(rho @ rho)
            population[i, j] = np.trace(n.full() @ rho)
            fidelity_to_squeezed[i, j], squeezing[i, j] = find_fidelity_and_squeezing(rho, rho_i)
        np.save('Figure 4/purity', purity)
        np.save('Figure 4/population', population)
        np.save('Figure 4/squeezing', squeezing)
        np.save('Figure 4/fidelity_to_squeezed', fidelity_to_squeezed)


def plot_populations():
    sigma_e = np.linspace(0.3, 10, 50)
    N = 50
    purity = np.load('Figure 4/purity.npy')
    squeezing = np.load('Figure 4/squeezing.npy')
    fidelity_to_squeezed = np.load('Figure 4/fidelity_to_squeezed.npy')
    polulation = np.load('Figure 4/population.npy')
    z = np.linspace(0, 49, 50, dtype=int)
    purity = purity[:, z != 26]
    squeezing = squeezing[:, z != 26]
    fidelity_to_squeezed = fidelity_to_squeezed[:, z != 26]
    sigma_e = sigma_e[z != 27]
    colors = ['blue', 'red', 'black']
    linestyles = ['-.', '-', '--']
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 2.5))
    for i in range(3):
        ax1.plot(sigma_e, purity[i, :], linewidth=3, color=colors[i], linestyle=linestyles[i])
        ax2.plot(sigma_e, squeezing[i, :], linewidth=3, color=colors[i], linestyle=linestyles[i])
        ax3.plot(sigma_e, fidelity_to_squeezed[i, :], linewidth=3, color=colors[i], linestyle=linestyles[i])

        ax1.set_xlabel('')
        # ax1.set_xscale('log')

        #
    for ax in (ax1, ax2, ax3):
        ax.set_xticks([])
        ax.set_xlim([0, 10])
        ax.set_yticks([])
        ax.set_xlim([0, 10])
    ax1.set_ylim([0.5, 1.02])
    ax2.set_ylim([0., 1.3])
    ax3.set_ylim([0.8, 1.005])
    plt.subplots_adjust(wspace=0.2, hspace=0.2)

    ax1.legend(['         ', '         ', '         '], fontsize=15)
    plt.savefig('populations.svg')
    plt.show()


def plt_wgnr(rho, f, axes):
    qutip.plot_wigner(qutip.Qobj(rho), fig=f, ax=axes, cmap='bwr',colorbar=True)
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xlabel('')
    axes.set_ylabel('')
    axes.set_xlim([-3, 3])
    axes.set_ylim([-5, 5])
    axes.set_title('')
    # axes.set_aspect('equal')


def calculate_fig_4_data():
    N_se = 50
    rho_u = qutip.ket2dm(qutip.basis(N, 0))
    for i, w in enumerate(np.linspace(0.3, 10, N_se)):
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))
        A = quad(epsilon, 0, T)[0]
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5
        rho = calculate_final_density_matrix(epsilon, rho_u)
        np.save('density_mats/Vacuum/eps_index=' + str(i), rho)
        print(i)
    rho_u = qutip.ket2dm(qutip.basis(N, 1))
    for i, w in enumerate(np.linspace(0.3, 10, N_se)):
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))
        A = quad(epsilon, 0, T)[0]
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5
        rho = calculate_final_density_matrix(epsilon, rho_u)
        np.save('density_mats/Fock/eps_index=' + str(i), rho)
        print(i)

    psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
    psi_i = psi_i / psi_i.norm()
    rho_u = qutip.ket2dm(psi_i)
    for i, w in enumerate(np.linspace(0.3, 10, N_se)):
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2))
        A = quad(epsilon, 0, T)[0]
        epsilon = lambda t: np.exp(-(t - t0 - 1) ** 2 / (2 * w ** 2)) / A * 1.5 * 5
        rho = calculate_final_density_matrix(epsilon, rho_u)
        np.save('density_mats/Cat/eps_index=' + str(i), rho)
        print(i)
    times_5_gain()


def modes_vs_width_epsilon_vacuum(times, widths, u, rho_u, s, num_modes=10):
    dt = times[1] - times[0]
    u = [u(t) for t in times]
    occupations = np.zeros([len(widths), num_modes])
    intensities = np.zeros([len(widths), len(times)])
    for i, w in enumerate(widths):
        epsilon = lambda t: np.exp(-(t - t0) ** 2 / (2 * w ** 2))
        A = quad(epsilon, 0, T)[0]
        epsilon = lambda t: np.exp(-(t - t0) ** 2 / (2 * w ** 2)) / A
        g1, _ = g1_t(times, epsilon, u, rho_u)
        print(i)
        occupations[i, :] = get_occupation_of_modes(g1 * dt, num_modes)

        np.save('maps/pulsed_drive/occupation_' + str(s), occupations)
        np.save('maps/pulsed_drive/intensities_' + str(s), intensities)
    return occupations


def modes_vs_width_epsilon(times, widths, u, rho_u, s, num_modes=2):
    dt = times[1] - times[0]
    u = [u(t) for t in times]
    occupations = np.zeros([len(widths), num_modes])
    intensities = np.zeros([len(widths), len(times)])
    for i, w in enumerate(widths):
        epsilon = lambda t: np.exp(-(t - t0) ** 2 / (2 * w ** 2))
        A = quad(epsilon, 0, T)[0]
        epsilon = lambda t: np.exp(-(t - t0) ** 2 / (2 * w ** 2)) / A
        g1, g1_input_dependant = g1_t(times, epsilon, u, rho_u)

        input_mode_1, input_mode_2 = get_most_occupied_mode_function(times, g1_input_dependant * dt, num_modes=2)
        print(conj(input_mode_1.T) @ g1 @ input_mode_1)
        occupations[i, 0] = (conj(input_mode_1.T) @ g1 @ input_mode_1 * dt ** 2)[0][0]
        occupations[i, 1] = (conj(input_mode_2.T) @ g1 @ input_mode_2 * dt ** 2)[0][0]
        print(occupations[i, :])

        np.save('maps/pulsed_drive/occupation_' + str(s), occupations)
        np.save('maps/pulsed_drive/intensities_' + str(s), intensities)
    return occupations


def calculate_final_density_matrix(epsilon, rho_u):
    times = np.linspace(0, T, 300)
    dt = times[1] - times[0]
    sigma_u = 1
    N = np.size(rho_u.full()[0, :])
    number = qutip.num(N)
    u_pulse = lambda t: np.exp(-(t - t0) ** 2 / (2 * sigma_u ** 2)) / (np.pi * sigma_u ** 2) ** .25
    u = [u_pulse(t) for t in times]
    # plt.plot([epsilon(t) for t in times])
    # plt.plot([u_pulse(t) for t in times])
    #
    # plt.show()
    g1, g1_input_dependant = g1_t(times, epsilon, u, rho_u)
    if qutip.expect(number, rho_u) == 0:
        print('vacuum')
        v = get_most_occupied_mode_function(times, g1 * dt)

    else:
        v = get_most_occupied_mode_function(times, g1 * dt)
    v_array = np.array([v(t) for t in times])
    print('occupation without Bloch messiah')
    print(conj(v_array.T) @ g1 @ v_array * dt ** 2)
    times = np.linspace(0, T, 1000)
    u = np.array([u_pulse(t) for t in times])

    dt = times[2] - times[1]
    v_arr = np.array([v(t) for t in times])
    np.save('v_func', v_arr)

    # plt.plot(times, u)
    # plt.plot(times, v_arr)
    # plt.show()
    coef_array = find_coefficiants_for_transformation_arrays(times, u, v, epsilon)
    rho_array = [rho_u, ket2dm(basis(N, 0)), ket2dm(basis(N, 0))]
    rho_f = transform_quantum_state(rho_array, coef_array)
    print('with BMT')
    print(qutip.expect(qutip.num(N), rho_f))
    return rho_f


def plot_figure_3():
    pad = 0.01
    # 3a
    f, axes = plt.subplots(2, 3, figsize=(6, 4))
    widths = np.logspace(-1, 1, 50)
    delays = np.linspace(-3, 3, 50)
    occupations_1 = np.load('Figure 3/3a/occupation_first_mode.npy')
    occupations_2 = np.load('Figure 3/3a/occupation_second_mode.npy')
    axes[0, 0].pcolormesh(widths, delays, occupations_1.T, cmap='inferno', shading='gouraud')
    axes[1, 0].pcolormesh(widths, delays, occupations_1.T / (occupations_1.T + occupations_2.T), cmap='inferno',
                          shading='gouraud')
    axes[1, 0].set_xscale('log')
    axes[0, 0].set_xscale('log')
    m0 = axes[0, 0].collections[0]
    m1 = axes[1, 0].collections[0]
    m1.set_clim(0.6, 1)
    cb1 = plt.colorbar(m0, ax=axes[0, 0], pad=pad)
    cb2 = plt.colorbar(m1, ax=axes[1, 0], pad=pad)
    # cb1.ax.set_yticks([])
    # cb2.ax.set_yticks([])

    # Figure 3b
    widths = np.logspace(-1, 1, 31)
    detunings = np.linspace(-2, 2, 31)
    occupations_1 = np.load('Figure 3/3b/occupation_first_mode.npy')
    occupations_2 = np.load('Figure 3/3b/occupation_second_mode.npy')
    axes[0, 1].pcolormesh(widths, detunings, occupations_1.T, cmap='inferno', shading='gouraud')
    axes[1, 1].pcolormesh(widths, detunings, occupations_1.T / (occupations_1.T + occupations_2.T), cmap='inferno',
                          shading='gouraud')
    axes[1, 1].set_xscale('log')
    axes[0, 1].set_xscale('log')
    m0 = axes[0, 1].collections[0]
    m1 = axes[1, 1].collections[0]
    m0.set_clim(2, 3)
    m1.set_clim(0.6, 0.8)
    cb1 = plt.colorbar(m0, ax=axes[0, 1], pad=pad)
    cb2 = plt.colorbar(m1, ax=axes[1, 1], pad=pad)
    # cb1.ax.set_yticks([])
    # cb2.ax.set_yticks([])

    with open('Figure 3/3c/TWPA_Deltas_and_xis.pkl', 'rb') as f:
        d = pickle.load(f)
        N = d['N']
        Deltas = d['Deltas']
        mode1_occupation = d['mode1_occupation']
        mode2_occupation = d['mode2_occupation']
        xis = d['xis']

        norm = mcolors.LogNorm(vmin=mode1_occupation.min(), vmax=mode1_occupation.max())
        axes[0, 2].pcolormesh(xis, Deltas, mode1_occupation, cmap='inferno', shading='gouraud', norm=norm)
        axes[1, 2].pcolormesh(xis, Deltas, mode1_occupation / (mode1_occupation + mode2_occupation),
                              cmap='inferno',
                              shading='gouraud')

        m0 = axes[0, 2].collections[0]
        m1 = axes[1, 2].collections[0]
        m0.set_clim([0, 1])
        cb1 = plt.colorbar(m0, ax=axes[0, 2], pad=pad)
        cb2 = plt.colorbar(m1, ax=axes[1, 2], pad=pad)
        # cb1.ax.set_yticks([])
        # cb2.ax.set_yticks([])
    for i in range(2):
        for j in range(3):
            axes[i, j].set_xticks([])
            axes[i, j].set_yticks([])

    # plt.savefig('fig_3.pdf')
    plt.show()


def calc_g1_omega(F, G, omegas, u, rho_u):
    Nw = len(omegas)
    u = np.array([u(omega) for omega in omegas])
    N = int(np.size(rho_u) ** 0.5)
    a = destroy(N)
    dw = omegas[1] - omegas[0]
    corr_n = qutip.expect(a.dag() * a, rho_u)
    corr_aa = qutip.expect(a * a, rho_u)
    corr_a_dag_a_dag = qutip.expect(a.dag() * a.dag(), rho_u)
    F_u_array = np.array([sum([F[i, k] * u[k] * dw for k in range(Nw)]) for i in range(Nw)])
    G_u_array = np.array([sum([G[i, k] * u[k] * dw for k in range(Nw)]) for i in range(Nw)])

    g1_omega = np.zeros([Nw, Nw], dtype=complex)
    g1_omega_input_dependant = np.zeros([Nw, Nw], dtype=complex)
    for i, w1 in enumerate(omegas):
        for j, w2 in enumerate(omegas):
            g1_omega[i, j] = corr_n * conj(F_u_array[i]) * F_u_array[j] + \
                             corr_a_dag_a_dag * conj(F_u_array[i]) * conj(G_u_array[j]) + \
                             corr_aa * F_u_array[j] * G_u_array[i] + corr_n * G_u_array[i] * conj(G_u_array[j]) + \
                             sum([G[i, k] * conj(G[j, k]) for k in range(Nw)])
            g1_omega_input_dependant[i, j] = corr_n * conj(F_u_array[i]) * F_u_array[j] + \
                                             corr_a_dag_a_dag * conj(F_u_array[i]) * conj(G_u_array[j]) + \
                                             corr_aa * F_u_array[j] * G_u_array[i] + corr_n * G_u_array[i] * conj(
                G_u_array[j])
    # f,(ax1,ax2) = plt.subplots(1,2)
    #
    # ax1.pcolormesh(omegas, omegas, np.abs(g1_omega), cmap='inferno')
    # ax2.pcolormesh(omegas, omegas, np.abs(g1_omega_input_dependant), cmap='inferno')
    #
    # plt.show()
    return g1_omega, g1_omega_input_dependant


def calculate_figure_3_data(a=False, b=False, c=False):
    '''

    :return:
    '''
    if a:
        # Figure 3a
        sigma_u = 1
        u = lambda t: np.exp(-(t - t0) ** 2 / (2 * sigma_u ** 2)) / (np.pi * sigma_u ** 2) ** .25
        times = np.linspace(0, T, 200)
        widths = np.logspace(-1, 1, 50)
        delays = np.linspace(-3, 3, 50)
        widths = [1]
        delays = [-1]
        rho_u = qutip.ket2dm(basis(N_se, 1))
        dt = times[1] - times[0]
        u = [u(t) for t in times]
        occupations_1 = np.zeros([len(widths), len(delays)])
        occupations_2 = np.zeros([len(widths), len(delays)])

        for i, w in enumerate(widths):
            print('i=' + str(i))
            for j, delta_t in enumerate(delays):
                print('j=' + str(j))
                epsilon = lambda t: np.exp(-(t - t0 - delta_t) ** 2 / (2 * w ** 2))
                A = quad(epsilon, 0, T)[0]
                print(A)
                epsilon = lambda t: np.exp(-(t - t0 - delta_t) ** 2 / (2 * w ** 2)) / A
                g1, g1_input_dependant = g1_t(times, epsilon, u, rho_u)
                input_mode_1, input_mode_2 = get_most_occupied_mode_function(times, g1_input_dependant * dt,
                                                                             num_modes=2)
                # plt.plot(times,input_mode_1)
                # plt.plot(times,input_mode_2)
                # plt.show()
                np.save('mode_1', input_mode_1)
                occupations_1[i, j] = (conj(input_mode_1.T) @ g1 @ input_mode_1 * dt ** 2)[0][0]
                occupations_2[i, j] = (conj(input_mode_2.T) @ g1 @ input_mode_2 * dt ** 2)[0][0]
                print(occupations_1[i, j])
                print(occupations_2[i, j])
            # np.save('Figure 3/3a/occupation_first_mode', occupations_1)
            # np.save('Figure 3/3a/occupation_second_mode', occupations_2)

    # Figure 3b
    if b:
        rho_u = qutip.ket2dm(basis(N_se, 1))

        # It includes the pump power and the nonlinearity of the medium.
        coupling = 20
        # Inverse group velocity of the pump beam
        A = 1
        # Inverse group velocity of the A beam
        B = 1
        # Inverse group velocity of the C beam
        C = 1

        # z_start gives the beginning of the crystal and z_stop its end
        # (z_start should remain at 0 since the analytic formulas in the code
        # assume a waveguide starting at 0.)
        z_start = 0
        z_stop = 2

        ###########################
        ## Evaluation parameters ##
        ###########################

        # Wavelength range to consider where the process is centered at 0. The values
        # have to cover all excited frequencies. 
        # (In some analysis scripts it is assumed  that w_start and w_stop are
        # symmetric about zero.)
        w_start = -5
        w_stop = 5
        # Sampling steps for the frequency degree of freedom
        w_steps = 100

        # Sampling steps for the propagation over the length of the crystal
        z_steps = 500
        # Width of the Gaussian pump field amplitude (not intensity) in sigma. 
        pump_widths = np.logspace(-1, 1, 31)
        # Directory to save the data in
        save_directory = "Data/FGs"
        detunings = np.linspace(-2, 2, 31)
        occupations_1 = np.zeros([len(pump_widths), len(detunings)])
        occupations_2 = np.zeros([len(pump_widths), len(detunings)])
        for i, pump_width in enumerate(pump_widths):
            print('i=' + str(i))

            for j, detuning in enumerate(detunings):
                print('j=' + str(j))
                sigma_u = 1
                u = lambda omega: np.exp(-(omega - detuning) ** 2 / (2 * sigma_u ** 2)) / (np.pi * sigma_u ** 2) ** .25

                F, G = calc_PDC(coupling, w_start, w_stop, w_steps, z_start, z_stop, z_steps, A, B, C, pump_width,
                                save_directory)
                omegas = np.linspace(w_start, w_stop, w_steps)
                dw = omegas[1] - omegas[0]

                g1, g1_input_dependant = calc_g1_omega(F, conj(G), omegas, u, rho_u)
                input_mode_1, input_mode_2 = get_most_occupied_mode_function(omegas, g1_input_dependant * dw,
                                                                             num_modes=2)
                # plt.plot(omegas,input_mode_1)
                # plt.plot(omegas,input_mode_2)
                # plt.show()

                occupations_1[i, j] = (conj(input_mode_1.T) @ g1 @ input_mode_1 * dw ** 2)[0][0]
                occupations_2[i, j] = (conj(input_mode_2.T) @ g1 @ input_mode_2 * dw ** 2)[0][0]
                print(occupations_1[i, j])
                print(occupations_2[i, j])
            np.save('Figure 3/3b/occupation_first_mode', occupations_1)
            np.save('Figure 3/3b/occupation_second_mode', occupations_2)


def calculate_figure_2_data():
    sigma_u = 1
    u = lambda t: np.exp(-(t - t0) ** 2 / (2 * sigma_u ** 2)) / (np.pi * sigma_u ** 2) ** .25
    times = np.linspace(0, T, 200)
    widths = np.logspace(-1, 1, 30)
    rho_u = qutip.ket2dm(qutip.basis(N, 0))
    modes_vs_width_epsilon_vacuum(times, widths, u, rho_u, 'Vacuum')

    rho_u = qutip.ket2dm(qutip.basis(N, 1))
    modes_vs_width_epsilon(times, widths, u, rho_u, 'Fock1')

    rho_u = qutip.ket2dm(qutip.basis(N, 4))
    modes_vs_width_epsilon(times, widths, u, rho_u, 'Fock4')

    rho_u = qutip.ket2dm(qutip.coherent(N, 2))
    modes_vs_width_epsilon(times, widths, u, rho_u, 'coherent2')

    psi_i = (displace(N, 2) + displace(N, -2)) * basis(N, 0)
    psi_i = psi_i / psi_i.norm()
    rho_u = qutip.ket2dm(psi_i)
    modes_vs_width_epsilon(times, widths, u, rho_u, 'cat2')


def pictures():
    times = np.linspace(0, T, 1000)
    sigma = 1
    mode_1 = np.load('mode_1.npy')
    # times = np.linspace(-5, 5, 10000)
    f, ax1 = plt.subplots(1, 1, figsize=(7, 5))
    # ax1.plot(np.exp(-times ** 2 / (2 * sigma ** 2)) * sin(times / 2), color='#0000CD', linewidth=1)
    ax1.plot(mode_1, color='#0000CD', linewidth=1)

    ax1.set_xticks([])
    ax1.set_yticks([])
    plt.axis('off')
    plt.xlim([40,160])

    plt.savefig('pump.svg')
    plt.show()

def plot_figure_2():
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(7, 3))
    widths = np.logspace(-1, 1, 30)
    occupations_vacuum = np.load('maps/pulsed_drive/occupation_Vacuum.npy')

    for i in range(len(occupations_vacuum[1, :])):
        ax1.plot(widths, occupations_vacuum[:, i], linewidth=2)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.set_ylim(0.001, 2)
        ax1.set_xlim([0.1, 10])

    occupation_coherent2 = np.load('maps/pulsed_drive/occupation_coherent2.npy')
    occupation_cat2 = np.load('maps/pulsed_drive/occupation_cat2.npy')
    occupations_fock4 = np.load('maps/pulsed_drive/occupation_Fock4.npy')
    occupations_fock1 = np.load('maps/pulsed_drive/occupation_Fock1.npy')

    for i in range(2):
        ax2.plot(widths, occupations_fock1[:, i], linewidth=2, color='red')
        ax2.plot(widths, occupations_fock4[:, i], linewidth=2, color='blue')
        if i == 0:
            ax2.plot(widths, occupation_coherent2[:, i], linewidth=2, color='orange')
            ax2.plot(widths, occupation_cat2[:, i], color='black', linewidth=2, linestyle='dashed')

        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.legend(['$\\vert n=1\\rangle$', '$\\vert n=4\\rangle$', '$\\vert \\alpha=2\\rangle$',
                    '$\\vert \\alpha=2\\rangle + \\vert \\alpha=-2\\rangle $'], frameon=False)

    ax2.set_ylim(0.02, 10)
    ax2.set_xlim([0.1, 10])
    plt.savefig('temp_fig2.svg')

    plt.show()


if __name__ == '__main__':
    # pictures()
    # plot_figure_2()
    # calculate_figure_2_data()
    # calculate_figure_3_data(a=True)
    # plot_figure_3()
    # plot_wigner_grid()
    # plot_populations()
    # calculate_figure_4_lineplots()
    #
    # plot_populations()
    plot_figure_4()
    # times_5_gain()
    # plot_figure_4()

    # calculate_fig_4_data()
