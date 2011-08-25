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

	t = copy(t)
	for i in range(len(t)):
		if t[i] < begin:
			t[i] = 0
		else:
			if (t[i] >= begin) and (t[i] < end):
				t[i] = (sin((pi/2.) * (t[i] - begin)/(end - begin)))**2
			else:
				t[i] = 1
	return t


def build_waveform(waveform_description, awg_samplingrate = 24):
	"""
	waveform = [
	{"amplitude": 1, "length": 3.2, "rise": 0.2, "frequency": 3.4, "phase":0, "fall": 0}
	]
	"""

	waveform = []
	for i in waveform_description:
		length = i["length"]
		N = round(length * awg_samplingrate)
		rise = i["rise"]
		fall = i["fall"]
		amplitude = i["amplitude"]
		frequency = i["frequency"]
		phase = i["phase"]
		t = arange(0, length, 1/float(awg_samplingrate))
		env = smooth(t, 0, rise) * (1 - smooth(t, length - fall, length)) 
		tmp_waveform = -512 + floor(1024. * env * amplitude * sin(2 * pi * (phase/360. + frequency * t)))
		waveform = concatenate((waveform, tmp_waveform));
	return waveform;


def build_signal(waveform_description, awg_samplingrate = 24):
	"""
	waveform = [
	{"amplitude": 1, "length": 3.2, "rise": 0.2, "frequency": 3.4, "phase":0, "fall": 0}
	]
	"""

	waveform = []
	for i in waveform_description:
		length = i["length"]
		rise = i["rise"]
		fall = i["fall"]
		amplitude = i["amplitude"]
		frequency = i["frequency"]
		phase = i["phase"]

		t = arange(0, length, 1/float(awg_samplingrate))
		env = smooth(t, 0, rise) * (1 - smooth(t, length - fall, length)) 
		tmp_waveform = env * amplitude * sin(2 * pi * (phase/360. + frequency * t))
		waveform = concatenate((waveform, tmp_waveform));
	return waveform;


def E(m, k, h):
    return - h[0]*m[0] + - h[1]*m[1] + - h[2]*m[2] - k[0] * m[0]**2 - k[1] * m[1]**2 - k[2] * m[2]**2  - k[3] * (m[0]**2 * m[1]**2 + m[1]**2 * m[2]**2 + m[2]**2 * m[0]**2)



def LL(m, h):
	return - cross(m, h)


def G(m, h, a):
	return -a * cross(m, cross(m, h))


def close_selected(arr):
	for i in arr:
		close(i)

def plot_curves(num, l, size):
	for i in l:
		openGraph(i);
		arr = getp(ci(), 'array')
		close()
		n = sqrt(len(arr))
		arr = arr.reshape((n, n))
		x = linspace(size[0], size[1], n)
		y = linspace(size[2], size[3], n)
		contourgraph(x, y, arr, 1);
		p = getp(getp(cA(), 'children')[2], 'paths')
		p = p[0].to_polygons()[0]
		a, b = hsplit(p, 2)
		close()
		figure(num);
		plot(a, b);


def adjust_size(a, new_size, new_value = 0):
	return concatenate((a, linspace(new_value, new_value, new_size - len(a))))

def astroid():
	THETA, PHI = meshgrid(linspace(0, pi/2, 100), linspace(0, 2*pi, 100)) 
	R = 1/(fabs(cos(THETA))**(2/3.) + fabs(sin(THETA)) ** (2/3.))**(3/2.)
	X = R * sin(THETA) * cos(PHI)
	Y = R * sin(THETA) * sin(PHI)
	Z = R * cos(THETA)
	emm.mesh(X, Y, Z)


def process_contour():
	arr = getp(ci(), 'array')
	print "arr size", len(arr), sqrt(len(arr))
	arr = arr.reshape((sqrt(len(arr)), sqrt(len(arr))))
	contourgraph(arr, 1)
	p = getp(getp(cA(), 'children')[2], 'paths')
	p = p[0].to_polygons()[0]
	a, b = hsplit(p, 2)
	plot(a, b)
	rr = sqrt(a**2 + b**2)
	rr = rr/rr.max()
	th = arctan2(b, a)
	print "th size", len(th)
	graph(rr*cos(th), rr*sin(th))
	th = pi/2. - th
	plot(rr*cos(th), rr*sin(th))
	r = meshgrid(rr, rr)[0]
	theta, phi = meshgrid(th, linspace(0, 2*pi, len(th))) 
	x = r * sin(theta) * cos(phi)
	y = r * sin(theta) * sin(phi)
	z = r * cos(theta)
	emm.mesh(x, y, z)


def process_surface(num_fig, angle):
	figure(num_fig)
	curves = []
	maxval = 0
	for i in cls():
		curves.append([getp(i, 'xdata').tolist(), getp(i, 'ydata').tolist()])
		if len(curves[-1][0]) > maxval:
			maxval = len(curves[-1][0])

	for i in curves:
		i[0] = adjust_size(i[0], maxval)
		i[1] = adjust_size(i[1], maxval)

	rr = []
	th = []
	for i in curves:
		rr.append(sqrt(i[0]**2 + i[1]**2))
		rr[-1] = rr[-1]/rr[-1].max()
		th.append(arctan2(i[1], i[0]))
		th[-1] = pi/2. - th[-1]
		# graph(rr[-1]*cos(th[-1]), rr[-1]*sin(th[-1]))

	theta, phi = meshgrid(linspace(0, pi/2., len(th[0])), angle)
	x = rr * sin(th) * cos(phi)
	y = rr * sin(th) * sin(phi)
	z = rr * cos(th)
	emm.mesh(x, y, z)


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




