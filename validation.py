
import subprocess
import compareFPDC

minipill = 1

d75 = 1

if minipill:
    dif = 0
    nmin=10
    nmax=40
    for i in range(nmin,nmax+1):
        command = f'./fcpw-libigl-fpdc-example ../minipill.stl 0 200004 {i}'
        proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        proc.wait()
        command2 = f'./fcpw-libigl-fpdc-example ../minipill.stl 2 200004 {i}'
        proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        proc2.wait()
        dif+=compareFPDC.compare(1,1,1)
        print(f'\n{i} : dif = {dif}\n')
    print("Malha = Oblong")
    print(f"Diferença acumulada para res=10,11,12,...,40 = {dif}")

if d75:
    dif = 0
    nmin=10
    nmax=40
    for i in range(nmin,nmax+1):
        command = f'./fcpw-libigl-fpdc-example ../d75.stl 0 200004 {i}'
        proc = subprocess.Popen(command, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        proc.wait()
        command2 = f'./fcpw-libigl-fpdc-example ../d75.stl 2 200004 {i}'
        proc2 = subprocess.Popen(command2, shell=False, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        proc2.wait()
        dif+=compareFPDC.compare(1,0,1)
        print(f'\n{i} : dif = {dif}\n')

    print(dif)

