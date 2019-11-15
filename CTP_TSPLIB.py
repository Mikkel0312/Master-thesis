from ConstructionCTPTSPLIB import *
from ImproveCTP import ImproveCTP
from MDMC import MDMC


def tspfile_to_df(filename):

    with open(filename) as file:
        dimension = None
        start = None
        end = None
        lines = file.readlines()

        i = 0
        while not end:
            line = lines[i]
            if line.startswith('DIMENSION :'):
                dimension = int(line.split()[2])
            if line.startswith('NODE_COORD_SECTION'):
                start = i
            if line.startswith('EOF'):
                end = i
            i = i+1


        file.seek(0)
        df = pd.read_csv(
            file,
            skiprows=start+1,
            sep=' ',
            names=['number', 'x', 'y'],
            dtype={'number': str, 'x': np.float64, 'y': np.float64},
            nrows=dimension
        )

        if (df.iloc[-1,0] == 'EOF'):
            df = df[:-1]
        df = df.drop('number', axis =1)



        return df

filename = "kroB150.tsp"
df = tspfile_to_df(filename)
T = 1

V = 25
W = len(df)-V

rl_ca, r_ca, cpu = solve_CTP(df, T, V , W, reduce = False, plot = True)




ImproveCTP(df, T, V, r_ca, rl_ca, plot = True)

df_W = df[V:]
df_W.index = df_W.index.astype(int)

