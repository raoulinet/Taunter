def from_mercator(X, Y):
    return [-tanh(Y), sqrt(1 -tanh(Y)**2)*cos(X), sqrt(1 -tanh(Y)**2)*sin(X)]



def to_mercator(x, y, z):
    """
    Trip in macrospin planet
    """
    # X, Y = [], []
    # for i in range(len(z)):
    #     X.append(-math.atan2(z[i], y[i]))
    #    Y.append(-math.atanh(x[i]))
    return [- arctan2(z, y), - arctanh(x)]



def tsv2rec(fname):
    """
    Overload of csv2rec function
    """
    return csv2rec(fname, delimiter=" ")


def smooth(t, begin, end):
    if t < begin:
        return 0
    else:
        if (t >= begin) and (t < end):
            return (sin((pi/2.) * (t - begin)/(end - begin)))**2
        else:
            return 1


def build_waveform(waveform_description, awg_samplingrate = 24):
    """
      waveform = [
      {"amplitude": 1, "length": 3.2, "rise": 0.2, "frequency": 3.4, "phase":0},
      {"amplitude": 0, "length": 1, "rise": 0.0, "frequency": 0, "phase": 0},
      {"amplitude": 1, "length": 20, "rise": 0.2, "frequency": ,"phase": 0}
      ]
    """
    
    awg_samplingrate = 1.*awg_samplingrate
    waveform = []
    for i in range(len(waveform_description)):
        length = waveform_description[i]["length"]
        N = round(length * awg_samplingrate)
        rise = waveform_description[i]["rise"]
        amplitude = waveform_description[i]["amplitude"]
        frequency = waveform_description[i]["frequency"]
        phase = waveform_description[i]["phase"]
        for j in range(N):
            t = j/awg_samplingrate
            # env = smooth(t, 0, rise) * (1 - smooth(t, length - rise, length)) 
            waveform.append(512 + 511 * amplitude * sin(2 * pi * (phase/360. + frequency * t)))
    waveform[len(waveform) - 1] = 511
    return waveform;



def E(m, k, h):
    return 1 - 0.5*k[0]*m[0]**2 - 0.5*k[1]*m[1]**2 - 0.5*k[2]*m[2]**2 + h[0]*m[0] + h[1]*m[1] + h[2]*m[2]


def powow(x, n):
    _resultat = []
    for _index in x:
        _resultat.append(pow(_index, n))
    return _resultat



def loadlast(folder_name):
    """
    Load last tsv file
    in the folder "folder_name"
    """
    ldir = os.listdir(folder_name)
    ldir.sort()
    fname = ldir.pop()	
    print("Load " + fname)
    ldir = []
    return tsv2rec(fname) 



def differential(x, t):
    """
    Raw derivation
    """
    dx_dt = array(x).copy()
    for i in range(len(dx_dt) - 2):
        dx_dt[i] = (x[i+1] - x[i])/(t[i+1] - t[i])
    dx_dt[-2] = dx_dt[-3]
    dx_dt[-1] = dx_dt[-3]
    return dx_dt



def unfold(phi, modulo):
    """
    Reverse modulo
    """
    unfold = []
    j = 0
    last = 0
    for i in range(len(phi) - 1):
        j = round((phi[i] - last)/(modulo))
        last = phi[i] - j*modulo
        phi[i] = last
    return phi



def mvtxt():
    """
    Move last axes text to clic point
    """

    _ax = gca()
    _place = ginput()[0]
    _ax.texts[-1].set_position(_place)
    draw()



def fit_leastsq(func = (lambda coef, a: coef[0] + coef[1] * a), data_X = [], data_Y = [], init_coef = []):
    """
    Least square fit routine
    func = lambda coef, a: coef[0] + coef[1] * a # Target function
    """
    
    data_X = array(data_X)
    data_Y = array(data_Y)
    init_coef = array(init_coef)
    errfunc = lambda coef, X, Y: func(coef, X) - Y # Distance to the target function
    return scipy.optimize.leastsq(errfunc, init_coef[:], args=(data_X, data_Y))




