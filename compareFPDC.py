import pandas as pd
import numpy as np

def compare(i,cb,iter=0):
    divergence = 1
    if i:
        b = ''


        df = pd.read_table(b+"I2_fpdc.dmat", sep = '\s+',skiprows=1,header=None, index_col=False)
        df2 = pd.read_table(b+"I0_fpdc.dmat", sep = '\s+',skiprows=1,header=None, index_col=False)

        df4 = df!=df2
        divergence = 1 - (df4.sum() / df4.shape[0]).to_numpy()[0]
        print(f"Contagem de valores diferentes de -1: {(df.replace(-1,np.nan).count())} ")
        if iter:
            return divergence
        print("Diferença FPDC/FPDC2")
        print("Size = ",df4.shape[0])
        
        divergence = df4[df4 == True].count().to_numpy()[0]
        formatted_divergence = divergence/df4.count().to_numpy()[0]
        print(divergence)
        print(f"Divergencia {formatted_divergence}")
        print(f"Contagem de valores diferentes de -1: {(df.replace(-1,np.nan).count())} ")
        input("Enter para continuar ......")
        print(df,df2,df4)
        df['key'] = df.index
        df2['key'] = df2.index
        df5 = pd.merge(df,df2,on='key')
        df5 = df5[['key','0_x','0_y']]
        print(df5)
        df_diff = df != df2
        # Filter rows where there are differences
        df_diffrows= df[df_diff.any(axis=1)]
        print(df_diffrows) 
        print(df5[df_diff.any(axis=1)])
        df5[df_diff.any(axis=1)].to_csv('diffROWS.txt', sep = ' ', index = False, header=False)
        

    if cb:
        b = ''
        
        df = pd.read_table(b+'id_mode0cellb.txt', sep = " ", header=None)
        df = df.iloc[:,1:]
        df.columns = ['x','y','z','w']
        print('df: cb')
        print(df)

        df2 = df.groupby(['x','y','z'])

        # print('filled cells')
        # print(len(df2.groups))
        # a = df.max()

        # print('total cells')
        # print(a['x']*a['y']*a['z'])

        # print('ids on first cell')
        # a = (df2.nth[0,0]['x'].iloc[0],df2.nth[0,0]['y'].iloc[0],df2.nth[0,0]['z'].iloc[0])
        # print (f'first cell: {a}')
        # print(df2.groups[(a[0],a[1],a[2])])


        # print('number of ids per cell')
        # print(df2.size())


        df1 = pd.read_table(b+'id_mode1cellb.txt', sep = " ", header=None)
        df1 = df1.iloc[:,1:]
        df1.columns = ['x','y','z','w']

        df12 = df1.groupby(['x','y','z'])

        # print('filled cells')
        # print(len(df12.groups))
        # a = df1.max()

        # print('total cells')
        # print(a['x']*a['y']*a['z'])

        # print('ids on first cell')
        # a = (df12.nth[0,0]['x'].iloc[0],df12.nth[0,0]['y'].iloc[0],df12.nth[0,0]['z'].iloc[0])
        # print (f'first cell: {a}')
        # print(df12.groups[(a[0],a[1],a[2])])

        # print('number of ids per cell')

        # df13 = df12 == df2
        # # for name_of_group, contents_of_group in df12.group_keys(0,1,5):
        # #     print(name_of_group)
        # #     print(list(contents_of_group))
        # # #     input("")
        # # print('\n\n DF12')
        # # group = df12.get_group((0, 1, 5))
        # # print(group['w'].sort_values().to_list())
        # # print('\n\n DF2')
        # # group = df2.get_group((0, 1, 5))
        # # print(group['w'].sort_values().to_list())
        print(list(df12.groups.keys())[-1])
        print(list(df2.groups.keys())[-1])

        print(df.shape, df1.shape)
        return divergence

compare(1,0)