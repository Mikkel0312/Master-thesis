from ConstructionCTP import solve_CTP, generate_points
from ImproveCTP import ImproveCTP
from MDMC import MDMC
import pandas as pd
import math


def table_CA(V_list, W_list, file1=None, file2=None):
    list1 = []
    T = 1
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            m = V + W
            df_CA = generate_points(m)
            for c in range(3):
                for s in range(3):
                    if c == 0 and s == 0:
                        continue
                    routelen_CA, route_CA, cpu = solve_CTP(df_CA, T, V, W, method="CA", cust=c, stemline=s)
                    list1.append([V, W, T, c, s, len(route_CA), routelen_CA])

    df_CA_results = pd.DataFrame(list1, columns=["|V|", "|W|", "|T|", "d(j)", "D(j)", "nodes", "obj"])

    d_0_1 = 0
    d_0_2 = 0
    d_1_0 = 0
    d_1_1 = 0
    d_1_2 = 0
    d_2_0 = 0
    d_2_1 = 0
    d_2_2 = 0

    for idx, i in df_CA_results.iterrows():
        if df_CA_results.loc[idx, "d(j)"] == 0 and df_CA_results.loc[idx, "D(j)"] == 1:
            d_0_1 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 0 and df_CA_results.loc[idx, "D(j)"] == 2:
            d_0_2 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 1 and df_CA_results.loc[idx, "D(j)"] == 0:
            d_1_0 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 1 and df_CA_results.loc[idx, "D(j)"] == 1:
            d_1_1 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 1 and df_CA_results.loc[idx, "D(j)"] == 2:
            d_1_2 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 2 and df_CA_results.loc[idx, "D(j)"] == 0:
            d_2_0 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 2 and df_CA_results.loc[idx, "D(j)"] == 1:
            d_2_1 += df_CA_results.loc[idx, "obj"]
        if df_CA_results.loc[idx, "d(j)"] == 2 and df_CA_results.loc[idx, "D(j)"] == 2:
            d_2_2 += df_CA_results.loc[idx, "obj"]

    list2 = [d_0_1, d_0_2, d_1_0, d_1_1, d_1_2, d_2_0, d_2_1, d_2_2]
    list2_str = ["d_0_1", "d_0_2", "d_1_0", "d_1_1", "d_1_2", "d_2_0", "d_2_1", "d_2_2"]
    list3 = [i / 6 for i in list2]
    df_compare = pd.DataFrame(list3, index=list2_str, columns=["obj"])
    if file1 != None:
      print(df_compare.to_latex(index = False), file = open(file1, 'w'))
    if file2 != None:
      print(df_CA_results.to_latex(index=False), file = open(file2, 'w'))
    print(df_compare.to_latex(index = False))
    print(df_CA_results.to_latex(index=False))



def table_OFFR(V_list, W_list, k_list, iterations = 3, file1=None, file2 = None):
    T = 1
    list1 = []
    for V in V_list:
        W_list1 = [V*w for w in W_list]
        for W in W_list1:
            for k in k_list:
                rl_OFFR, cpu1 = 0,0
                for i in range(iterations):
                    m = V + W
                    df_CA = generate_points(m)
                    routelen_OFFR, route_1, cpu, a, b = solve_CTP(df_CA, T, V, W, method="OFFR", k=k, returnparameter=True)
                    rl_OFFR += routelen_OFFR
                    cpu1 += cpu
                rl_OFFR = rl_OFFR/iterations
                cpu1 = cpu1 / iterations
                list1.append([V, W, k, a, b, rl_OFFR, cpu1])

    OFFR_results = pd.DataFrame(list1, columns=["V", "W", "k", "a", "b", "obj", "CPU"])
    avg_k_list = []
    for k in k_list:
        k_reults = OFFR_results['k'] == k
        avg_k_list.append([k, OFFR_results.loc[k_reults, "a"].mean(), OFFR_results.loc[k_reults, "b"].mean(), OFFR_results.loc[k_reults, "obj"].mean(), OFFR_results.loc[k_reults, "CPU"].mean()])

    avg_k_df = pd.DataFrame(avg_k_list, columns = ["k", "a", "b","obj", "CPU"])

    print(OFFR_results.to_latex(index=False))
    print(avg_k_df.to_latex(index = False))

    if file1 != None:
        print(OFFR_results.to_latex(index=False), file=open(file1, 'w'))
    if file2 != None:
        print(avg_k_df.to_latex(index=False), file=open(file2, 'w'))



def table_ONR(V_list, W_list, k_list, iterations = 5, cust = 0, file1=None, file2 = None):
    T = 1
    list1 = []
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            for k in k_list:
                rl_ONR, cpu1 = 0, 0
                for i in range(iterations):
                    m = V + W
                    df_CA = generate_points(m)
                    routelen_ONR, route_1, cpu, a, b = solve_CTP(df_CA, T, V, W, method="ONR", k=k,cust = cust,
                                                                  returnparameter=True)
                    rl_ONR += routelen_ONR
                    cpu1 += cpu
                rl_ONR = rl_ONR / iterations
                cpu1 = cpu1 / iterations
                list1.append([V, W, k, a, b, rl_ONR, cpu1])

    ONR_results = pd.DataFrame(list1, columns=["V", "W", "k", "a", "b", "obj", "CPU"])
    avg_k_list = []
    for k in k_list:
        k_reults = ONR_results['k'] == k
        avg_k_list.append([k, ONR_results.loc[k_reults, "a"].mean(), ONR_results.loc[k_reults, "b"].mean(),
                           ONR_results.loc[k_reults, "obj"].mean(), ONR_results.loc[k_reults, "CPU"].mean()])

    avg_k_df = pd.DataFrame(avg_k_list, columns=["k", "a", "b", "obj", "CPU"])

    print(ONR_results.to_latex(index=False))
    print(avg_k_df.to_latex(index=False))

    if file1 != None:
        print(ONR_results.to_latex(index=False), file=open(file1, 'w'))
    if file2 != None:
        print(avg_k_df.to_latex(index=False), file=open(file2, 'w'))


def table_compare_regression(V_list, W_list, iterations = 2, file1=None):
    T = 1
    list1 = []
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            list_OFFR = []
            list_ONR = []
            for i in range(iterations):
                m = V + W
                df = generate_points(m)
                routelen_ONR, route_1, cpu_ONR, a_ONR, b_ONR = solve_CTP(df, T, V, W, method="ONR", k=20, cust=1,
                                                             returnparameter=True)
                routelen_OFFR, route_1, cpu_OFFR, a_OFFR, b_OFFR = solve_CTP(df, T, V, W, method="OFFR", k=5, cust=1,
                                                             returnparameter=True)
                list_OFFR.append([routelen_OFFR, cpu_OFFR, a_OFFR, b_OFFR])
                list_ONR.append([routelen_ONR, cpu_ONR, a_ONR, b_ONR])


            rl_OFFR = sum([i[0] for i in list_OFFR])/len(list_OFFR)
            cpu_OFFR = sum([i[1] for i in list_OFFR])/len(list_OFFR)
            a_OFFR = sum([i[2] for i in list_OFFR])/len(list_OFFR)
            b_OFFR = sum([i[3] for i in list_OFFR])/len(list_OFFR)

            rl_ONR = sum([i[0] for i in list_ONR]) / len(list_ONR)
            cpu_ONR = sum([i[1] for i in list_ONR]) / len(list_ONR)
            a_ONR = sum([i[2] for i in list_ONR]) / len(list_ONR)
            b_ONR = sum([i[3] for i in list_ONR]) / len(list_ONR)

            list1.append([V, W, a_OFFR, b_OFFR, rl_OFFR, cpu_OFFR, a_ONR, b_ONR, rl_ONR, cpu_ONR])



    results = pd.DataFrame(list1, columns=["V", "W", "a_OFFR", "b_OFFR", "obj_OFFR", "CPU_OFFR", "a_ONR", "b_ONR",
                                            "obj_ONR", "CPU_ONR"])

    list1.append(["mean", "", results['a_OFFR'].mean(), results['b_OFFR'].mean(), results['obj_OFFR'].mean(),
                  results['CPU_OFFR'].mean(), results['a_ONR'].mean(), results['b_ONR'].mean(), results['obj_ONR'].mean(), results['CPU_ONR'].mean()])

    results = pd.DataFrame(list1, columns=["V", "W", "a_OFFR", "b_OFFR", "obj_OFFR", "CPU_OFFR", "a_ONR", "b_ONR",
                                           "obj_ONR", "CPU_ONR"])
    if file1 != None:
        print(results.to_latex(index=False), file=open(file1, 'w'))
    print(results.to_latex(index=False))


def table_compare_RLA(V_list, W_list, iterations = 5, file1=None):
    T = 1
    list1 = []
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            sum_1, sum_2, sum_3, sum_4, sum_5 = [0, 0, 0, 0, 0]
            for i in range(iterations):
                m = V + W
                df = generate_points(m)
                r_Bea, ro, cpu_bea = solve_CTP(df, T, V, W, method="RLA:1")
                r_Dag, route_Dag, cpu_dag = solve_CTP(df, T, V, W, method="RLA:2")
                r_Chi, route_1, cpu_Chien = solve_CTP(df, T, V, W, method="RLA:3")
                r_Kwo_1, route_1, cpu_Kwon1 = solve_CTP(df, T, V, W, method="RLA:4")
                r_Kwo_2, route_1, cpu_Kwon2 = solve_CTP(df, T, V, W, method="RLA:5")
                sum_1 += r_Bea
                sum_2 += r_Dag
                sum_3 += r_Chi
                sum_4 += r_Kwo_1
                sum_5 += r_Kwo_2
            avg_Bea = sum_1 / iterations
            avg_Dag = sum_2 / iterations
            avg_Chi = sum_3 / iterations
            avg_Kwo_1 = sum_4 / iterations
            avg_Kwo_2 = sum_5 / iterations
            list1.append([V, W, avg_Bea, avg_Dag, avg_Chi, avg_Kwo_1, avg_Kwo_2])
    results = pd.DataFrame(list1, columns=["V", "W", "avg_Bea", "avg_Dag", "avg_Chi", "avg_Kwo_1", "avg_Kwo_2"])

    list1.append(["mean", "", results['avg_Bea'].mean(), results['avg_Dag'].mean(), results['avg_Chi'].mean(), results['avg_Kwo_1'].mean(), results['avg_Kwo_2'].mean(), ])
    results = pd.DataFrame(list1, columns=["V", "W", "avg_Bea", "avg_Dag", "avg_Chi", "avg_Kwo_1", "avg_Kwo_2"])

    print(results.to_latex(index=False))

    if file1 != None:
        print(results.to_latex(index=False), file=open(file1, 'w'))




def table_parameters_ImproveCTP(V_list, W_list, iterations = 3, file1=None, file2 = None):
    T = 1
    list1 = []
    list2 = []
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            m = V + W
            df = generate_points(m)

            routelen_CA, route_CA, cpu_CA = solve_CTP(df, T, V, W, method="CA", cust=0, stemline=1)
            improved_CA_len_1, route, cpu_1 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=1)
            improved_CA_len_2, route, cpu_2 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=2)
            improved_CA_len_3, route, cpu_3 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=3)
            improved_CA_len_n, route, cpu_4 = ImproveCTP(df, T, V, route_CA, routelen_CA)
            [routelen_CA, improved_CA_len_1, cpu_1, improved_CA_len_2, cpu_2, improved_CA_len_3, cpu_3, improved_CA_len_n,
            cpu_4] = [round(i, 2) for i in [routelen_CA, improved_CA_len_1, cpu_1, improved_CA_len_2, cpu_2, improved_CA_len_3, cpu_3, improved_CA_len_n,
                 cpu_4]]
            list1.append([V, W, routelen_CA, improved_CA_len_1, cpu_1, improved_CA_len_2, cpu_2, improved_CA_len_3, cpu_3, improved_CA_len_n,
                 cpu_4])
    for V in V_list:
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            m = V + W
            df = generate_points(m)
            routelen_CA, route_CA, cpu_CA = solve_CTP(df, T, V, W, method="CA", cust=0, stemline=1)

            improved_CA_len_25, route, cpu_25 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=1, maxiters=math.ceil(0.25 * len(route_CA)))
            improved_CA_len_50, route, cpu_50 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=1,
                                                          maxiters=math.ceil(0.50 * len(route_CA)))
            improved_CA_len_75, route, cpu_75 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=1,
                                                          maxiters=math.ceil(0.75 * len(route_CA)))
            improved_CA_len_100, route, cpu_100 = ImproveCTP(df, T, V, route_CA, routelen_CA, n=1,
                                                          maxiters=math.ceil(len(route_CA)))
            [V, W, routelen_CA, improved_CA_len_25, improved_CA_len_50, improved_CA_len_75, improved_CA_len_100] = [round(i,2) for i in [V, W, routelen_CA, improved_CA_len_25, improved_CA_len_50, improved_CA_len_75, improved_CA_len_100]]
            list2.append([V, W, routelen_CA, improved_CA_len_25, improved_CA_len_50, improved_CA_len_75, improved_CA_len_100])

    results_alpha = pd.DataFrame(list1, columns = ["|V|", "|W|" , "CA_obj", "obj_1", "cpu_1" , "obj_2" , "cpu_2", "obj_3", "cpu_3", "obj_4", "cpu_4"])


    list1.append(["mean", "", results_alpha["CA_obj"].mean(), results_alpha["obj_1"].mean(), results_alpha["cpu_1"].mean(), results_alpha["obj_2"].mean(), results_alpha["cpu_2"].mean(), results_alpha["obj_3"].mean(), results_alpha["cpu_3"].mean(), results_alpha["obj_4"].mean(), results_alpha["cpu_4"].mean()])
    results_alpha = pd.DataFrame(list1,
                                 columns=["|V|", "|W|", "CA_obj", "obj_1", "cpu_1", "obj_2", "cpu_2", "obj_3", "cpu_3",
                                          "obj_4", "cpu_4"])
    print(results_alpha.to_latex(index = False))

    results_maxiter = pd.DataFrame(list2, columns = ["|V|" , "|W|", "obj_CA",  "0.25", "0.50", "0.75", "1"])

    list2.append(["mean", " ", results_maxiter["obj_CA"].mean(),results_maxiter["0.25"].mean(), results_maxiter["0.50"].mean(), results_maxiter["0.75"].mean(), results_maxiter["1"].mean()])
    results_maxiter = pd.DataFrame(list2, columns=["|V|" , "|W|", "obj_CA",  "0.25", "0.50", "0.75", "1"])

    print(results_maxiter.to_latex(index = False))

    if file1 != None:
        print(results_alpha.to_latex(index=False), file=open(file1, 'w'))
    if file2 != None:
        print(results_maxiter.to_latex(index=False), file=open(file2, 'w'))


def table_compare_methods(T, V_list, W_list, iterations=5, MDMC_iter=5, smallscale=True, file1=None):
    list_results = []

    for v in V_list:
        V = v

        W_list1 = [V * w for w in W_list]
        for W in W_list1:

            list_ca, list_reg, list_rla, list_MDMC = [], [], [], []
            for i in range(iterations):
                m = W + V
                df = generate_points(m)
                rl_ca, r_ca, cpu_ca = solve_CTP(df, T, V, W, method="CA", cust=0, stemline=1)
                if not smallscale:
                    rl_reg, r_reg, cpu_reg = solve_CTP(df, T, V, W, method="OFFR", k=25, cust =1)
                    list_reg.append([rl_reg, cpu_reg])
                rl_rla, r_rla, cpu_rla = solve_CTP(df, T, V, W, method="RLA:3")
                list_ca.append([rl_ca, cpu_ca])
                list_rla.append([rl_rla, cpu_rla])

                rl_MDMC, cpu_MDMC = 0, 0
                best_MDMC = 99999999
                lol_MDMC = []
                for j in range(MDMC_iter):
                    rl_MDMC_cur, r_MDMC, cpu_MDMC_cur = MDMC(df, T, V)
                    if rl_MDMC_cur < best_MDMC:
                        best_MDMC = rl_MDMC_cur
                    rl_MDMC += rl_MDMC_cur
                    cpu_MDMC += cpu_MDMC_cur
                    lol_MDMC.append([rl_MDMC_cur, cpu_MDMC_cur])
                cpu_MDMC = cpu_MDMC / MDMC_iter
                rl_MDMC = rl_MDMC / MDMC_iter
                list_MDMC.append([rl_MDMC, cpu_MDMC, best_MDMC])
                print(V, V_list, W, W_list1, i)
            avg_rl_ca = round(sum([i[0] for i in list_ca]) / iterations,2)
            avg_cpu_ca = round(sum([i[1] for i in list_ca]) / iterations, 2)
            if not smallscale:
                avg_rl_reg = round(sum([i[0] for i in list_reg]) / iterations,2)
                avg_cpu_reg = round(sum([i[1] for i in list_reg]) / iterations, 2)

            avg_rl_rla = round(sum([i[0] for i in list_rla]) / iterations,2)
            avg_cpu_rla = round(sum([i[1] for i in list_rla]) / iterations, 2)

            avg_rl_MDMC = round(sum([i[0] for i in list_MDMC]) / iterations,2)
            avg_cpu_MDMC = round(sum([i[1] for i in list_MDMC]) / iterations, 2)

            if not smallscale:
                best = min(avg_rl_MDMC, avg_rl_ca, avg_rl_reg, avg_rl_rla)
            else:
                best = min(avg_rl_MDMC, avg_rl_ca, avg_rl_rla)
            gap_ca = round((avg_rl_ca - best) / avg_rl_ca*100, 2)
            if not smallscale:
                gap_reg = round((avg_rl_reg - best) / avg_rl_reg*100, 2)
            gap_rla = round((avg_rl_rla - best) / avg_rl_rla*100, 2)
            gap_MDMC = round((avg_rl_MDMC - best) / avg_rl_MDMC*100, 2)
            if not smallscale:
                list_results.append(
                    [V, W, avg_rl_ca, avg_cpu_ca, gap_ca, avg_rl_reg, avg_cpu_reg, gap_reg, avg_rl_rla,
                     avg_cpu_rla, gap_rla, avg_rl_MDMC, avg_cpu_MDMC, gap_MDMC])
            else:
                list_results.append(
                    [V, W, avg_rl_ca, avg_cpu_ca, gap_ca, avg_rl_rla, avg_cpu_rla, gap_rla, avg_rl_MDMC,
                     avg_cpu_MDMC, gap_MDMC])
    if not smallscale:
        df_results = pd.DataFrame(list_results,
                                  columns=["V", "W", "avg_rl_ca", "avg_cpu_ca", "gap_ca", "avg_rl_reg",
                                           "avg_cpu_reg", "gap_reg", "avg_rl_rla", "avg_cpu_rla", "gap_rla",
                                           "avg_rl_MDMC", "avg_cpu_MDMC", "gap_MDMC"])
        avg_1 = df_results["avg_rl_ca"].mean()
        avg_2 = df_results["avg_cpu_ca"].mean()
        avg_3 = df_results["gap_ca"].mean()
        avg_4 = df_results["avg_rl_reg"].mean()
        avg_5 = df_results["avg_cpu_reg"].mean()
        avg_6 = df_results["gap_reg"].mean()
        avg_7 = df_results["avg_rl_rla"].mean()
        avg_8 = df_results["avg_cpu_rla"].mean()
        avg_9 = df_results["gap_rla"].mean()
        avg_10 = df_results["avg_rl_MDMC"].mean()
        avg_11 = df_results["avg_cpu_MDMC"].mean()
        avg_12 = df_results["gap_MDMC"].mean()
        [avg_1, avg_2, avg_3, avg_4, avg_5, avg_6, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12] = [round(i,2) for i in [avg_1, avg_2, avg_3, avg_4, avg_5, avg_6, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12]]
        list_results.append(
            [ " ", " ", avg_1, avg_2, avg_3, avg_4, avg_5, avg_6, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12])

        df_results = pd.DataFrame(list_results,
                                  columns=["V", "W", "rl_ca", "cpu_ca", "gap_ca", "el_reg",
                                           "cpu_reg", "gap_reg", "avg_rl_rla", "cpu_rla", "gap_rla",
                                           "rl_MDMC", "cpu_MDMC", "gap_MDMC"])

    else:
        df_results = pd.DataFrame(list_results,
                                  columns=["V", "W", "avg_rl_ca", "avg_cpu_ca", "gap_ca", "avg_rl_rla",
                                           "avg_cpu_rla", "gap_rla",
                                           "avg_rl_MDMC", "avg_cpu_MDMC", "gap_MDMC"])
        avg_1 = df_results["avg_rl_ca"].mean()
        avg_2 = df_results["avg_cpu_ca"].mean()
        avg_3 = df_results["gap_ca"].mean()
        avg_7 = df_results["avg_rl_rla"].mean()
        avg_8 = df_results["avg_cpu_rla"].mean()
        avg_9 = df_results["gap_rla"].mean()
        avg_10 = df_results["avg_rl_MDMC"].mean()
        avg_11 = df_results["avg_cpu_MDMC"].mean()
        avg_12 = df_results["gap_MDMC"].mean()
        [avg_1, avg_2, avg_3, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12] = [round(i, 2) for i in [avg_1, avg_2, avg_3, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12]]
        list_results.append([" ", " ", avg_1, avg_2, avg_3, avg_7, avg_8, avg_9, avg_10, avg_11, avg_12])

        df_results = pd.DataFrame(list_results,
                              columns=["V", "W", "avg_rl_ca", "avg_cpu_ca", "gap_ca", "avg_rl_rla",
                                       "avg_cpu_rla", "gap_rla",
                                       "avg_rl_MDMC", "avg_cpu_MDMC", "gap_MDMC"])


    print(df_results.to_latex(index = False))
    if file1 != None:
        print(df_results.to_latex(index = False), file=open(file1, 'w'))




def table_compare_methods_improve(T, V_list, W_list, iterations=5, MDMC_iter=5, file1=None, file2 = None):
    list_results1 = []
    list_results2 = []
    for v in V_list:
        V = v
        W_list1 = [V * w for w in W_list]
        for W in W_list1:
            list_ca, list_reg, list_rla, list_MDMC = [], [], [], []
            for i in range(iterations):
                m = W + V
                df = generate_points(m)
                rl_ca, r_ca, cpu_ca = solve_CTP(df, T, V, W, method="CA", cust=0, stemline=1)
                irl_ca, ir_ca, icpu_ca = ImproveCTP(df, T, V, r_ca, rl_ca, n=1)
                rl_reg, r_reg, cpu_reg = solve_CTP(df, T, V, W, method="ONR", k=50)
                irl_reg, ir_reg, icpu_reg = ImproveCTP(df, T, V, r_reg, rl_reg, n=1)
                rl_rla, r_rla, cpu_rla = solve_CTP(df, T, V, W, method="RLA:2")
                irl_rla, ir_rla, icpu_rla = ImproveCTP(df, T, V, r_reg, rl_reg, n=1)

                list_ca.append([irl_ca, cpu_ca+icpu_ca, 100 *(rl_ca-irl_ca)/rl_ca])
                list_reg.append([irl_reg, cpu_reg + icpu_reg, 100 *(rl_reg-irl_reg)/rl_reg])
                list_rla.append([irl_rla, cpu_rla + icpu_rla, 100 * (rl_rla-irl_rla)/rl_rla])

                rl_MDMC, cpu_MDMC = 0, 0
                irl_MDMC, icpu_MDMC = 0,0
                best_MDMC = 99999999
                for j in range(MDMC_iter):
                    rl_MDMC_cur, r_MDMC, cpu_MDMC_cur = MDMC(df, T, V)
                    irl_MDMC_cur, ir_MDMC, icpu_MDMC_cur = ImproveCTP(df, T, V, r_MDMC, rl_MDMC_cur, n=1)

                    if irl_MDMC_cur < best_MDMC:
                        best_MDMC = irl_MDMC_cur
                    rl_MDMC += rl_MDMC_cur
                    cpu_MDMC += cpu_MDMC_cur+icpu_MDMC_cur
                    irl_MDMC += irl_MDMC_cur
                cpu_MDMC = cpu_MDMC / MDMC_iter
                rl_MDMC = rl_MDMC / MDMC_iter
                irl_MDMC = irl_MDMC / MDMC_iter
                list_MDMC.append([rl_MDMC, cpu_MDMC, 100 * (rl_MDMC - irl_MDMC)/rl_MDMC])
                print(V, V_list, W, W_list1, i)

            avg_rl_ca = round(sum([i[0] for i in list_ca]) / iterations,2)
            avg_cpu_ca = round(sum([i[1] for i in list_ca]) / iterations, 2)
            avg_improve_ca = round(sum([i[2] for i in list_ca]) / iterations, 2)

            avg_rl_reg = round(sum([i[0] for i in list_reg]) / iterations,2)
            avg_cpu_reg = round(sum([i[1] for i in list_reg]) / iterations, 2)
            avg_improve_reg = round(sum([i[2] for i in list_reg]) / iterations, 2)

            avg_rl_rla = round(sum([i[0] for i in list_rla]) / iterations,2)
            avg_cpu_rla = round(sum([i[1] for i in list_rla]) / iterations, 2)
            avg_improve_rla = round(sum([i[2] for i in list_rla]) / iterations, 2)

            avg_rl_MDMC = round(sum([i[0] for i in list_MDMC]) / iterations,2)
            avg_cpu_MDMC = round(sum([i[1] for i in list_MDMC]) / iterations, 2)
            avg_improve_MDMC = round(sum([i[2] for i in list_MDMC]) / iterations, 2)

            best = min(avg_rl_MDMC, avg_rl_ca, avg_rl_reg, avg_rl_rla)

            gap_ca = round(((avg_rl_ca - best) / avg_rl_ca)* 100 , 2)
            gap_reg = round(((avg_rl_reg - best) / avg_rl_reg)* 100, 2)
            gap_rla = round(((avg_rl_rla - best) / avg_rl_rla)* 100, 2)
            gap_MDMC = round(((avg_rl_MDMC - best) / avg_rl_MDMC)* 100, 2)

            list_results1.append(
                [T, V, W, avg_rl_ca, avg_cpu_ca, gap_ca, avg_improve_ca, avg_rl_reg, avg_cpu_reg, gap_reg, avg_improve_reg])
            list_results2.append([T, V, W, avg_rl_rla, avg_cpu_rla,
                 gap_rla, avg_improve_rla, avg_rl_MDMC, avg_cpu_MDMC, gap_MDMC, avg_improve_MDMC])

    df_results1 = pd.DataFrame(list_results1,
                              columns=["T", "V", "W", "avg_rl_ca", "avg_cpu_ca", "gap_ca", "avg_improve_ca", "avg_rl_reg", "avg_cpu_reg",
                                       "gap_reg", "avg_improve_reg"])
    df_results2 = pd.DataFrame(list_results2, columns = ["T", "V", "W", "avg_rl_rla", "avg_cpu_rla", "gap_rla", "avg_improve_rla", "avg_rl_MDMC", "avg_cpu_MDMC",
                                       "gap_MDMC", "avg_improve_MDMC"])

    avg_1 = df_results1["avg_rl_ca"].mean()
    avg_2 = df_results1["avg_cpu_ca"].mean()
    avg_3 = df_results1["gap_ca"].mean()
    avg_3i = df_results1["avg_improve_ca"].mean()
    avg_4 = df_results1["avg_rl_reg"].mean()
    avg_5 = df_results1["avg_cpu_reg"].mean()
    avg_6 = df_results1["gap_reg"].mean()
    avg_6i = df_results1["avg_improve_reg"].mean()
    avg_7 = df_results2["avg_rl_rla"].mean()
    avg_8 = df_results2["avg_cpu_rla"].mean()
    avg_9 = df_results2["gap_rla"].mean()
    avg_9i = df_results2["avg_improve_rla"].mean()
    avg_10 = df_results2["avg_rl_MDMC"].mean()
    avg_11 = df_results2["avg_cpu_MDMC"].mean()
    avg_12 = df_results2["gap_MDMC"].mean()
    avg_12i = df_results2["avg_improve_MDMC"].mean()

    list_results1.append(
        [" ", " ", " ", avg_1, avg_2, avg_3, avg_3i, avg_4, avg_5, avg_6, avg_6i])
    list_results2.append([" ", " ", " ",avg_7, avg_8, avg_9, avg_9i, avg_10, avg_11, avg_12, avg_12i])

    df_results1 = pd.DataFrame(list_results1,
                              columns=["T", "V", "W", "obj_ca", "cpu", "GAP", "improve", "obj_reg", "cpu", "gap", "improve"])
    df_results2 = pd.DataFrame(list_results2, columns = ["T", "V", "W", "avg_rl_rla", "avg_cpu_rla", "gap_rla", "avg_improve_rla", "avg_rl_MDMC", "avg_cpu_MDMC",
                                       "gap_MDMC", "avg_improve_MDMC"])

    print(df_results1.to_latex(index = False))
    print(df_results2.to_latex(index = False))
    if file1 != None:
        print(df_results1.to_latex(index = False), file=open(file1, 'w'))
    if file2 != None:
        print(df_results2.to_latex(index = False), file=open(file2, 'w'))


# results for the CA parameter settings.
# V_list = [20, 30, 40]
# W_list = [1,3]
# iterations = 3
# file1,file2 = "CA1.txt", "CA2.txt"
# table_CA(V_list, W_list, file1 = file1, file2 = file2)

# results for the OFFR parameter settings
# V_list = [30,40,50]
# W_list = [1,3]
# k_list = [10,25,50]
# file1 = "OFFR1.txt"
# file2 = "OFFR2.txt"
# table_OFFR(V_list, W_list, k_list, iterations = 2, cust=1, file1 = file1, file2 = file2)

#results for the ONR parameter settings
# V_list = [30,40,50]
# W_list = [1,3]
# k_list = [10,25,50, 100, 200]
# file1 = "ONR1.txt"
# file2 = "ONR2.txt"
# table_ONR(V_list, W_list, k_list, iterations = 2, cust=1, file1 = file1, file2 = file2)

# V_list = [25,50,75]
# W_list = [1,3]
# file1 = "regcompare.txt"
# table_compare_regression(V_list, W_list, iterations = 5, file1=file1)


#results for the RLA parameter settings
#V_list = [25,50,75]#,50]
#W_list = [1,3]
#file1 = "RLA.txt"
#table_compare_RLA(V_list, W_list, iterations = 2, file1=file1)

#results for the ImproveCTP parameter settings
# V_list = [25,50,75]
# W_list = [1, 3]
# file1 = "Imp1.txt"
# file2 = "Imp2.txt"
# table_parameters_ImproveCTP(V_list, W_list, file1=file1, file2 = file2)

#results for comparison in small scale
# T = 1
# V_list = [5,10,15]
# W_list = [1,2]
# table_compare_methods(T, V_list, W_list, iterations=5, MDMC_iter=5, smallscale=True, file1="CompareSmallscaleT1.txt")


#results for comparison in medium scale
#V_list = [20, 30, 40]
#W_list = [1, 2, 3]
#T = 1

#table_compare_methods(T, V_list, W_list, iterations=3, MDMC_iter=5, smallscale=True, file1="CompareMediumscaleT1.txt")


#results for comparison in large scale
#V_list = [25,50,75,100]
#W_list = [1,2,3,4]
#T = 1

#table_compare_methods(T, V_list, W_list, iterations=5, MDMC_iter=5, smallscale=False, file1="CompareLargesacleT1.txt")