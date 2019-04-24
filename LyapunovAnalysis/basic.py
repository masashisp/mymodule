
def CalcDky(lamb):
    print('test')
    Dky = 0.0
    d=len(lamb)
    k=0
    hp=0
    for p in range(d):
        hp += lamb[p]
        if hp <= 0.0:
            hp -= lamb[p]
            k=p
            break
        else:
            k=p+1
    print(k)
    if k==0:
        Dky=0.0
    elif k==d:
        Dky=d
    else:
       Dky = k + hp/abs(lamb[k])

    return Dky
