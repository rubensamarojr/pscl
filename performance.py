import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

d75 = 0

if d75:
    dif = 0
    nmin=10
    nmax=110
    nParts = 4
    print('Testing Pre Process...')
    bufferIt1 = np.zeros((11,2))
    fillIt1 = np.zeros((11,2))
    bufferMean = np.zeros((11,2))
    fillMean = np.zeros((11,2))
    totalppMean = np.zeros((11,2))
    for c,i in enumerate(range(nmin,nmax+1,10)):
        
        
        command = f'./fcpw-libigl-fpdc-example ../d75.stl 0 {nParts} {i}'

        proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
        proc.wait()
        command2 = f'./fcpw-libigl-fpdc-example ../d75.stl 1 {nParts} {i}'

        proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        proc2.wait()
        print(f"i={i}")     
    
    for c,i in enumerate(range(nmin,nmax+1,10)):
        for j in range(2):
            df = pd.read_csv(f'output/bm{j}nP{nParts}res{i}_fpdc.txt', sep= ';',header=None)
            
            bufferIt1[c,j] = df.iloc[2,1]
            fillIt1[c,j] = df.iloc[3,1]
            bufferMean[c,j] = df.iloc[6,1]
            fillMean[c,j] = df.iloc[7,1]
            totalppMean[c,j] = df.iloc[6,1] + df.iloc[7,1]

            print(f"c={c}, i={i}, j={j}")



    plt.plot(range(nmin,nmax+1,10),fillMean[:,0], label = "fillMode0", color = 'red')
    plt.plot(range(nmin,nmax+1,10),fillMean[:,1], label = "fillMode1", color = 'red', linestyle='dashed')

    plt.plot(range(nmin,nmax+1,10),bufferMean[:,0], label = "bufferMode0", color = 'blue')
    plt.plot(range(nmin,nmax+1,10),bufferMean[:,1], label = "bufferMode1", color = 'blue', linestyle='dashed')

    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,0], label = "totalMode0", color = 'orange')
    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,1], label = "totalMode1", color = 'orange', linestyle='dashed')
    plt.legend(ncols=2)
    plt.show()

d75_q = 0

if d75_q:

    dif = 0
    nmin = 10
    nmax = 60
    nParts = 1000004
    print('Testing Pre Process...')
    bufferIt1 = np.zeros((12,3))
    fillIt1 = np.zeros((12,3))
    bufferMean = np.zeros((12,3))
    fillMean = np.zeros((12,3))
    totalMean = np.zeros((12,3))
    for c,i in enumerate(range(nmin,nmax+1,5)):
        
        
       command = f'./fcpw-libigl-fpdc-example ../d75.stl 0 {nParts} {i}'

       proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
       proc.wait()
       command2 = f'./fcpw-libigl-fpdc-example ../d75.stl 1 {nParts} {i}'

       proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
       proc2.wait()

       
       print(f"c={c}, i={i}")
    #    command3 = f'../oldfpdc/fcpw-libigl-fpdc-example ../d75.stl {nParts} {i}'
    #    proc3 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    #    proc3.wait()
    #    print(f"i={i}")     
    
    for c,i in enumerate(range(nmin,nmax+1,5)):
        for j in range(2):
        #     if (j<2):  df = pd.read_csv(f'output/bm{j}nP{nParts}res{i}_fpdc.txt', sep= ';',header=None)
        #     if (j==2): df = pd.read_csv(f'output/bmnP{nParts}res{i}_oldfpdc.txt', sep= ';',header=None)
            df = pd.read_csv(f'output/bm{j}nP{nParts}res{i}_fpdc.txt', sep= ';',header=None)
            bufferMean[c,j] = df.iloc[-2,1]
            fillMean[c,j] = df.iloc[-1,1] 
            totalMean[c,j] = df.iloc[-2,1] 


    aa = list(range(nmin,nmax+1,5))

    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,2], label = "Query Flat FPDC", color = 'blue', linestyle='dotted')
    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,0], label = "Query Optimized Old FPDC", color = 'blue')
    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,1], label = "Query New FPDC", color = 'blue', linestyle='dashed')
    
    plt.plot(range(nmin,nmax+1,5),totalMean[:-1,2], label = "Total Flat FPDC", color = 'red', linestyle='dotted')
    plt.plot(range(nmin,nmax+1,5),totalMean[:-1,0], label = "Total Optimized Old FPDC", color = 'red')
    plt.plot(range(nmin,nmax+1,5),totalMean[:-1,1], label = "Total New FPDC", color = 'red', linestyle='dashed')


    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,0], label = "totalMode0", color = 'orange')
    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,1], label = "totalMode1", color = 'orange', linestyle='dashed')
    plt.legend(ncols=2)
    plt.show()

minipill_q = 1

if minipill_q:

    dif = 0
    nmin = 10
    nmax = 80
    nParts = 10000004
    print('Testing Pre Process...')
    bufferIt1 = np.zeros((16,3))
    fillIt1 = np.zeros((16,3))
    bufferMean = np.zeros((16,3))
    fillMean = np.zeros((16,3))
    totalppMean = np.zeros((16,3))
    for c,i in enumerate(range(nmin,nmax+1,5)):
        
        
       command = f'./fcpw-libigl-fpdc-example ../minipill.stl 0 {nParts} {i}'

       proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
       proc.wait()
       command2 = f'./fcpw-libigl-fpdc-example ../minipill.stl 1 {nParts} {i}'

       proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
       proc2.wait()
    #    command3 = f'./fcpw-libigl-fpdc-example ../minipill.stl 2 {nParts} {i}'
    #    proc3 = subprocess.Popen(command3, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    #    proc3.wait()
    #    print(f"i={i}")     
    
    for c,i in enumerate(range(nmin,nmax+1,5)):
        for j in range(3):
            df = pd.read_csv(f'output/bm{j}nP{nParts}res{i}_fpdc.txt', sep= ';',header=None)
          
            bufferMean[c,j] = df.iloc[-2,1]
            fillMean[c,j] = df.iloc[-1,1] 
            totalppMean[c,j] = df.iloc[6,1] + df.iloc[7,1] 

            print(f"c={c}, i={i}, j={j}")

    aa = list(range(nmin,nmax+1,5))

    # plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,2], label = "Query Flat FPDC", color = 'blue', linestyle='dotted')
    # plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,0], label = "Query Optimized Old FPDC", color = 'blue')
    # plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,1], label = "Query New FPDC", color = 'blue', linestyle='dashed')
    
    # plt.plot(range(nmin,nmax+1,5),fillMean[:-1,2], label = "Total Flat FPDC", color = 'red', linestyle='dotted')
    # plt.plot(range(nmin,nmax+1,5),fillMean[:-1,0], label = "Total Optimized Old FPDC", color = 'red')
    # plt.plot(range(nmin,nmax+1,5),fillMean[:-1,1], label = "Total New FPDC", color = 'red', linestyle='dashed')
    
    plt.plot(range(nmin,nmax+1,5),totalppMean[:-1,2], label = "Total Flat FPDC", color = 'red', linestyle='dotted')
    plt.plot(range(nmin,nmax+1,5),totalppMean[:-1,0], label = "Total Optimized Old FPDC", color = 'red')
    plt.plot(range(nmin,nmax+1,5),totalppMean[:-1,1], label = "Total New FPDC", color = 'red', linestyle='dashed')


    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,0], label = "totalMode0", color = 'orange')
    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,1], label = "totalMode1", color = 'orange', linestyle='dashed')
    plt.legend(ncols=2)
    plt.show()

minipill = 0
if minipill:
    dif = 0
    nmin=10
    nmax=160
    nParts = 4
    print('Testing Pre Process...')
    bufferIt1 = np.zeros((16,2))
    fillIt1 = np.zeros((16,2))
    bufferMean = np.zeros((16,2))
    fillMean = np.zeros((16,2))
    totalppMean = np.zeros((16,2))
    for c,i in enumerate(range(nmin,nmax+1,5)):
        
        command = f'./fcpw-libigl-fpdc-example ../minipill.stl 0 {nParts} {i}'

        proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
        proc.wait()
        command2 = f'./fcpw-libigl-fpdc-example ../minipill.stl 1 {nParts} {i}'

        proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
        proc2.wait()

        command3 = f'./fcpw-libigl-fpdc-example ../minipill.stl 2 {nParts} {i}'

        proc3 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
               
        proc3.wait()

        
        print(f"i={i}")            

    for c,i in enumerate(range(nmin,nmax+1,5)):
        
        for j in range(3):
            df = pd.read_csv(f'output/bm{j}nP{nParts}res{i}_fpdc.txt', sep= ';',header=None)
          
            bufferMean[c,j] = df.iloc[-2,1]
            fillMean[c,j] = df.iloc[-1,1] 

            print(f"c={c}, i={i}, j={j}")

        print(f"c={c}, i={i}, j={j}")




    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,2], label = "Query Flat FPDC", color = 'blue', linestyle='dotted')
    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,0], label = "Query Optimized Old FPDC", color = 'blue')
    plt.plot(range(nmin,nmax+1,5),bufferMean[:-1,1], label = "Query New FPDC", color = 'blue', linestyle='dashed')
    
    plt.plot(range(nmin,nmax+1,5),fillMean[:-1,2], label = "Total Flat FPDC", color = 'red', linestyle='dotted')
    plt.plot(range(nmin,nmax+1,5),fillMean[:-1,0], label = "Total Optimized Old FPDC", color = 'red')
    plt.plot(range(nmin,nmax+1,5),fillMean[:-1,1], label = "Total New FPDC", color = 'red', linestyle='dashed')
    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,0], label = "totalMode0", color = 'orange')
    # plt.plot(range(nmin,nmax+1,10),totalppMean[:,1], label = "totalMode1", color = 'orange', linestyle='dashed')
    plt.legend(ncols=2)
    plt.show()
