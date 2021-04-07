def generating_plots(Rp):

    fig, (plot1ax, plot2ax) = plt.subplots(1, 2)

    print("\n***Plotting desulpurization degree versus C_SO2_0***\n")

    start_time_alpha = timeit.default_timer()

    # x-axis : concentration of C_SO2_0
    # y-axis: desulphurization degree
    # keeping t_mean fixed
    # varying alpha

    alpha_list = [1, 2, 3, 4]
    t_mean_fixed = 4

    for alpha in alpha_list:

        conc_points = 50

        C_SO2_0 = np.linspace(200, 1400, conc_points)

        D = np.zeros(conc_points)

        for i in range(conc_points):
            X_CaO = find_X_CaO(alpha, ppmv_to_molm3(C_SO2_0[i], T1), t_mean_fixed, Rp)
            D[i] = desulphurization_degree(alpha, X_CaO)

        plot1ax.plot(C_SO2_0, D, label=alpha_label(alpha))

        print("Done with alpha = ", alpha, " after t = ", timeit.default_timer() - start_time_alpha, " seconds")



    print("\n***Plotting desulpurization degree versus t_mean***\n")

    start_time_t = timeit.default_timer()

    alpha_fixed = 2
    conc_list = [200, 500, 1000, 1400]

    for C_SO2_0 in conc_list:

        t_points = 50

        t = np.linspace(0.001, 8, t_points)

        D = np.zeros(t_points)

        for i in range(t_points):
            X_CaO = find_X_CaO(alpha_fixed, ppmv_to_molm3(C_SO2_0, T1), t[i], Rp)
            D[i] = desulphurization_degree(alpha_fixed, X_CaO)
        
        plot2ax.plot(t, D, label=C_SO2_label(C_SO2_0))
        
        print("Done with C_SO2 = ", C_SO2_0, " after t = ", timeit.default_timer() - start_time_t)




    
    

    plot1ax.set(xlabel="Inlet concentration of SO2 $C(SO_2)_0$ ppmv", ylabel=desulph_str)
    plot1ax.set_title("Keeping " + r"$\overline{t}_{mean}$ = " + str(t_mean_fixed))

    plot2ax.set(xlabel="Mean residence time " + r"$\overline{t}$ [h]", ylabel=desulph_str)
    plot2ax.set_title("Keeping " + r"$\alpha = $" + str(alpha_fixed))

    plot1ax.legend()
    plot2ax.legend()

    fig.suptitle("With diameter = " + str(Rp*2*10**3) + " [mm]")

    plt.show()