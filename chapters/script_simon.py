import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import time
import pandas as pd

import one_d_problem_physics as opp


def compare_K_w_script():
    # import moisture_transport_edit.K_w_calc
    # import moisture_transport_edit.P_suc_calc
    from moisture_transport_edit import P_suc_calc as P_suc_andi
    from moisture_transport_edit import K_w_calc as K_w_andi
    prob = opp.one_d_problem(res=20, sim_time=24, material='original', length=.1)

    n = prob.n
    A = prob.A
    free_saturation = prob.free_saturation
    pore_size = prob.pore_size

    print("andis K_w_calc(P_suc([1,10,100]))")
    print(K_w_andi(P_suc_andi(np.array([1,10,100]),free_saturation, pore_size), n, A, pore_size))

    print("my K_w_calc(P_suc([1,10,100]))")
    print(prob.K_w_calc(prob.P_suc_calc(np.array([1,10,100]))))


def benchmark_script():
    print('########## BENCHMARK SCRIPT #############')
    ress = [20, 40, 80, 160, 320]
    methods = ['BDF', 'LSODA']
    vectorized = [True, False]
    #mats = ['original', 'Ziegel1']
    mats = ['Ziegel1']

    results = {'resolution': ress}
    for m in methods:
        for v in vectorized:
            for mat in mats:
                timings = []
                for r in ress:
                    from scipy.integrate import solve_ivp
                    my_problem = opp.one_d_problem(res=r, sim_time=24, material=mat, length=.1)
                    w0 = my_problem.w[:]
                    t0 = 0  # Start time in hours
                    tf = my_problem.sim_time  # End time in hours
                    t_eval = np.linspace(t0, tf, 100)

                    print(f'Solving the differential equation... resolution: {r}, vectorized: {v}, method: {m}, material: {mat}')
                    start_time = time.time()
                    sol = solve_ivp(my_problem.dwdt_calc_vec, (t0, tf), w0, t_eval=t_eval, method=m, vectorized=v, dense_output=False, atol=1e-7,
                                    rtol=1e-5)
                    # sol = solve_ivp(my_problem.dwdt_calc_vec, (t0, tf), w0, t_eval=t_eval, method='Radau', vectorized=True, dense_output=False, atol=1e-7, rtol=1e-5)
                    del sol
                    t = time.time() - start_time
                    print(f'time elapsed: {t:.3f} s')
                    timings.append(t)

                key = m + '_' + mat + '_v' + str(v)
                results[key] = timings
    return results

def plot_benchmark(filename):
    df = pd.read_csv(filename, sep=';')
    df.set_index('resolution', inplace=True)
    df.plot(loglog=True, use_index=True)
    plt.ylabel('run time in s')
    plt.grid()
    plt.show()


def aufnahme_koef_WIP(sol, dx, index=20):
    return (sol.y[1:-1, index+1].sum() - sol.y[1:-1, index].sum()) * dx / (np.sqrt(sol.t[index+1]) - np.sqrt(sol.t[index]))


def validation_sctipt1():
    # imbibition...
    # --------------

    val_problem = opp.one_d_problem(res=160, sim_time=100, material='AAC_A4_mod', init_w=0, length=.1, w_west=1, w_east=0)
    val_problem.fluid_flow_west = True
    val_problem.fluid_flow_east = False
    val_problem.vapour_flow_west = False
    val_problem.vapour_flow_east = False

    val_problem.liquid_conduction = True
    val_problem.vapour_diffusion = True

    t0 = 0  # Start time in hours
    tf = val_problem.sim_time  # End time in hours
    t_eval = np.linspace(t0, tf, 100)
    w0 = val_problem.w[:]

    print('Solving the differential equation...')
    start_time = time.time()
    val_sol = solve_ivp(val_problem.dwdt_calc, (t0, tf), w0, t_eval=t_eval, method='Radau',
                    vectorized=False, dense_output=False, atol=1e-7, rtol=1e-5)
    print(f'time elapsed: {(time.time() - start_time):.3f} s')

    plt.plot(val_problem.x, val_sol.y[1:-1,1], label="1 h")
    plt.plot(val_problem.x, val_sol.y[1:-1,24], label="24 h")
    plt.plot(val_problem.x, val_sol.y[1:-1,48], label=">48 h")
    plt.legend()
    plt.grid()
    plt.xlabel('specimen length in $m$')
    plt.ylabel('water content in $kg/m^3$')
    plt.show()



def validation_sctipt2():
    # drying...
    # ------------

    # calculate P_suc and w for 23°C, 50%rH:
    P_suc_west = opp.p_suc_kelvin_calc(0.5)


    vval_problem = opp.one_d_problem(res=80, sim_time=100, material='AAC_A4_mod', init_w=0, length=.1, w_west=0, w_east=0)
    w_sat = vval_problem.free_saturation
    w_west = vval_problem.w_calc(P_suc_west) / w_sat
    val_problem = opp.one_d_problem(res=80, sim_time=1000, material='AAC_A4_mod_dry', init_w=310/w_sat, length=.1, w_west=w_west, w_east=0)
    # val_problem = opp.one_d_problem(res=30, sim_time=1000, material='AAC_A4_mod', init_w=310/350, length=.1, w_west=0, w_east=0)


    val_problem.fluid_flow_west = False
    val_problem.fluid_flow_east = False
    val_problem.vapour_flow_west = True
    val_problem.vapour_flow_east = False

    val_problem.liquid_conduction = True
    val_problem.vapour_diffusion = True

    t0 = 0  # Start time in hours
    tf = val_problem.sim_time  # End time in hours
    t_eval = np.linspace(t0, tf, 1000)
    w0 = val_problem.w[:]

    print('Solving the differential equation...')
    start_time = time.time()
    val_sol = solve_ivp(val_problem.dwdt_calc, (t0, tf), w0, t_eval=t_eval, method='Radau',
                    vectorized=False, dense_output=False, atol=1e-7, rtol=1e-5)
    print(f'time elapsed: {(time.time() - start_time):.3f} s')

    plt.plot(val_problem.x, val_sol.y[1:-1, 1], label="1 h")
    plt.plot(val_problem.x, val_sol.y[1:-1, 16], label="16 h")
    plt.plot(val_problem.x, val_sol.y[1:-1, 42], label="42 h")
    plt.plot(val_problem.x, val_sol.y[1:-1, 120], label="120 h")
    plt.plot(val_problem.x, val_sol.y[1:-1, 257], label="257 h")
    plt.plot(val_problem.x, val_sol.y[1:-1, 664], label="664 h")
    plt.legend()
    plt.grid()
    plt.xlabel('specimen length in $m$')
    plt.ylabel('water content in $kg/m^3$')
    plt.show()





######################################
#              MAIN
######################################
def main():
    validation_sctipt2()
    # plot_benchmark('benchmarks_RK_Radau.csv')
    benchmarks = False
    single_problem = False
    plotting = False

    # ----------------------
    #       Benchmarks
    # ---------------------
    if benchmarks:
        start_time = time.time()
        results = benchmark_script()
        print(f'total benchmark time: {(time.time() - start_time):.3f} s')
        pd.DataFrame(results).to_csv('benchmarks.csv', sep=';')


    # ---------------------------
    #       Solve Single Problem
    # ---------------------------
    if single_problem:
        #my_problem = opp.one_d_problem(res=40, sim_time=48, material='original', length=.1)
        #TODO: fix init_w and check pycharm todos
        my_problem = opp.one_d_problem(res=80, sim_time=100, material='AAC_A4_mod', init_w = 0, length=.1, w_west=1, w_east=0)

        my_problem.fluid_flow_west = True
        my_problem.fluid_flow_east = False
        my_problem.vapour_flow_west = False
        my_problem.vapour_flow_east = False

        my_problem.liquid_conduction = True
        my_problem.vapour_diffusion = True

        print(f'Material: {my_problem.material}')
        w0 = my_problem.w[:]

        t0 = 0  # Start time in hours
        tf = my_problem.sim_time  # End time in hours

        t_eval = np.linspace(t0, tf, 100)

        print('Solving the differential equation...')
        start_time = time.time()
        #sol = solve_ivp(my_problem.dwdt_calc, (t0, tf), w0, t_eval=t_eval, dense_output=False, atol=1e-7, rtol=1e-5)
        sol = solve_ivp(my_problem.dwdt_calc, (t0, tf), w0, t_eval=t_eval, method='Radau',
                        vectorized=False, dense_output=False, atol=1e-7, rtol=1e-5)
        print(f'time elapsed: {(time.time() - start_time):.3f} s')
        # print(f'w after one timestep (t = {sol.t[1]}):')
        # print(sol.y[:,1])
        # print(f'w final:')
        # print(sol.y[:,-1])

    # ----------------------
    #       Plots
    # ---------------------
    plt.rcParams['text.usetex'] = False
    if plotting:

        fig, axes = plt.subplots(ncols=3, figsize=(20,5))

        # Moisture content distribution in discretized specimen
        axes[0].grid()

        # Betrachtungsraum wählen
        axes[0].axis([0, 1.1*my_problem.length, 0, 1.1*my_problem.free_saturation])
        axes[0].set_xlabel('Discretized specimen in $m$')
        axes[0].set_ylabel('Moisture content for cell in $kg/m^3$')
        #axes[0].set_title(f'Moisture distribution in specimen')
        box_string = r'$\beta = $' + '(nA)\n' + r'$pore\_size = $' + f'{my_problem.pore_size}\n' + r'$n = $' + f'{my_problem.n}\n' + \
                     r'$A = $' + f'{my_problem.A}\n'
        #Adding text inside a rectangular box by using the keyword 'bbox'
        axes[0].annotate(box_string, xy=(0.02, 0.75), xycoords='axes fraction', bbox = dict(facecolor = 'orange', alpha = 0.25))

        # Moisture content in specimen over square root of time
        axes[1].grid()
        axes[1].set_xlabel('$\sqrt{Time}$ in $\sqrt{h}$')
        axes[1].set_ylabel('Total moisture content in $kg/m^2$')
        axes[1].axis([0, 1.1*np.sqrt(sol.t[-1]), 0, 1.1*max(sum(sol.y[:,0]), sum(sol.y[:,-1]))*my_problem.dx])

        # Relative moisture content in specimen over square root of time
        axes[2].grid()
        axes[2].set_xlabel('$\sqrt{Time}$ in $\sqrt{h}$')
        axes[2].set_ylabel('Relative moisture content in $\%$')
        #axes[2].set_title(f'Development of the relative moisture content in the specimen')
        axes[2].axis([0, 1.1*np.sqrt(sol.t[-1]), 0, 105])

        # two ways of making animated results
        animated_results(axes, my_problem, sol, anim_time = 3) # this one stays
        #animated_results_FuncAnimation(fig, axes, my_problem, sol, anim_time=3) # this one loops
        plt.show()


def animated_results(axes, problem, sol, anim_time=4):
    # animation parameters
    fps = 20
    total_frames = sol.y.shape[1]
    frame_freq = int(total_frames / (fps * anim_time)) + 1
    display_frames = [i * frame_freq for i in range(fps*anim_time) if i*frame_freq < sol.y.shape[1]]
    display_frames.append(sol.y.shape[1]-1)
    anim_pause = anim_time / (fps * sol.y.shape[1])

    # initial states
    line1, = axes[0].plot(problem.x, sol.y[1:-1,0])
    t_sqrt_axis = np.sqrt(sol.t)
    line2, = axes[1].plot([t_sqrt_axis[0]], [sum(sol.y[1:-1,0] * problem.dx)])
    line3, = axes[2].plot([t_sqrt_axis[0]], [sum(sol.y[1:-1,0] * problem.dx / (problem.free_saturation * problem.length) * 100)])

    # animation
    for i in display_frames:
        line1.set_ydata(sol.y[1:-1, i])
        line2.set_xdata(t_sqrt_axis[:i+1])
        line2.set_ydata([sum(sol.y[1:-1,j]) * problem.dx for j in range(i+1)])
        line3.set_xdata(t_sqrt_axis[:i+1])
        line3.set_ydata([sum(sol.y[1:-1,j]) * problem.dx / (problem.free_saturation * problem.length) * 100
                        for j in range(i + 1)])
        plt.pause(anim_pause)


def animated_results_FuncAnimation(fig, axes, problem, sol, anim_time=4):
    import matplotlib.animation as anim
    fps = 20
    total_frames = sol.y.shape[1]
    frame_freq = int(total_frames / (fps * anim_time)) + 1
    display_frames = [i * frame_freq for i in range(fps*anim_time) if i*frame_freq < sol.y.shape[1]]
    display_frames.append(sol.y.shape[1]-1)

    # initial states
    line1, = axes[0].plot(problem.x, sol.y[1:-1,0])
    t_sqrt_axis = np.sqrt(sol.t)
    line2, = axes[1].plot([t_sqrt_axis[0]], [sum(sol.y[1:-1,0] * problem.dx)])
    line3, = axes[2].plot([t_sqrt_axis[0]], [sum(sol.y[1:-1,0] * problem.dx / (problem.free_saturation * problem.length) * 100)])
    def my_animation(i):
        line1.set_ydata(sol.y[1:-1, i])
        line2.set_xdata(t_sqrt_axis[:i+1])
        line2.set_ydata([sum(sol.y[1:-1,j]) * problem.dx for j in range(i+1)])
        line3.set_xdata(t_sqrt_axis[:i+1])
        line3.set_ydata([sum(sol.y[1:-1,j]) * problem.dx / (problem.free_saturation * problem.length) * 100
                        for j in range(i + 1)])
        if i == sol.y.shape[1]-1:
            time.sleep(anim_time/2)

    ani = anim.FuncAnimation(fig, my_animation, frames=display_frames, interval=int(1000/fps))
    plt.show()




if __name__ == '__main__':
    pass
    main()

