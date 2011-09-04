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



def skim_array(arr, minima, maxima, col):

	to_skim = copy(arr)
	to_skim = to_skim.tolist()
	for i in range(len(arr)):
		for j in range(len(arr[i])):
			if arr[i][j] == maxima:
				to_skim[i][j] = (1, 0, 0, 0)
			else:
				if arr[i][j] == minima:
					to_skim[i][j] = (0, 1, 0, 0)
				else:
					to_skim[i][j] = col
	return to_skim



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


def mouse_fit():
	_points = ginput(2, show_clicks = False)
	return (_points[1][1] - _points[0][1])/(_points[1][0] - _points[0][0])


def fit_linear(data_X = [], data_Y = [], init_coef = []):
    """
    Least square fit routine
    func = lambda coef, a: coef[0] + coef[1] * a # Target function
    """
    func = (lambda coef, a: coef[0] + coef[1] * a) 
    data_X = array(data_X)
    data_Y = array(data_Y)
    init_coef = array(init_coef)
    errfunc = lambda coef, X, Y: func(coef, X) - Y # Distance to the target function
    return scipy.optimize.leastsq(errfunc, init_coef[:], args=(data_X, data_Y))


def fit_leastsq(func, data_X = [], data_Y = [], init_coef = []):
    """
    Least square fit routine
    func = lambda coef, a: coef[0] + coef[1] * a # Target function
    """
    
    data_X = array(data_X)
    data_Y = array(data_Y)
    init_coef = array(init_coef)
    errfunc = lambda coef, X, Y: func(coef, X) - Y # Distance to the target function
    return scipy.optimize.leastsq(errfunc, init_coef[:], args=(data_X, data_Y))




def griddata(x, y, z, binsize=0.01, retbin=True, retloc=True):
	"""
	Place unevenly spaced 2D data on a grid by 2D binning (nearest
	neighbor interpolation).

	Parameters
	----------
	x : ndarray (1D)
	The idependent data x-axis of the grid.
	y : ndarray (1D)
	The idependent data y-axis of the grid.
	z : ndarray (1D)
	The dependent data in the form z = f(x,y).
	binsize : scalar, optional
	The full width and height of each bin on the grid.  If each
	bin is a cube, then this is the x and y dimension.  This is
	the step in both directions, x and y. Defaults to 0.01.
	retbin : boolean, optional
	Function returns `bins` variable (see below for description)
	if set to True.  Defaults to True.
	retloc : boolean, optional
	Function returns `wherebins` variable (see below for description)
	if set to True.  Defaults to True.

	Returns
	-------
	grid : ndarray (2D)
	The evenly gridded data.  The value of each cell is the median
	value of the contents of the bin.
	bins : ndarray (2D)
	A grid the same shape as `grid`, except the value of each cell
	is the number of points in that bin.  Returns only if
	`retbin` is set to True.
	wherebin : list (2D)
	A 2D list the same shape as `grid` and `bins` where each cell
	contains the indicies of `z` which contain the values stored
	in the particular bin.

	Revisions
	---------
	2010-07-11  ccampo  Initial version
	"""
	# get extrema values.
	xmin, xmax = x.min(), x.max()
	ymin, ymax = y.min(), y.max()

	# make coordinate arrays.
	xi      = np.arange(xmin, xmax+binsize, binsize)
	yi      = np.arange(ymin, ymax+binsize, binsize)
	xi, yi = np.meshgrid(xi,yi)

	# make the grid.
	grid           = np.zeros(xi.shape, dtype=x.dtype)
	nrow, ncol = grid.shape
	if retbin: bins = np.copy(grid)

	# create list in same shape as grid to store indices
	if retloc:
		wherebin = np.copy(grid)
		wherebin = wherebin.tolist()

	# fill in the grid.
	for row in range(nrow):
		for col in range(ncol):
			xc = xi[row, col]    # x coordinate.
			yc = yi[row, col]    # y coordinate.

			# find the position that xc and yc correspond to.
			posx = np.abs(x - xc)
			posy = np.abs(y - yc)
			ibin = np.logical_and(posx < binsize/2., posy < binsize/2.)
			ind  = np.where(wbin == True)[0]

			# fill the bin.
			bin = z[ibin]
			if retloc: wherebin[row][col] = ind
			if retbin: bins[row, col] = bin.size
			if bin.size != 0:
				binval         = np.median(bin)
				grid[row, col] = binval
			else:
				grid[row, col] = np.nan   # fill empty bins with nans.

			# return the grid
	if retbin:
		if retloc:
			return grid, bins, wherebin
		else:
			return grid, bins
	else:
		if retloc:
			return grid, wherebin
		else:
			return grid



def cmap_map(function,cmap):
	""" Applies function (which should operate on vectors of shape 3:
	[r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.

	light_jet = cmap_map(lambda x: x/2+0.5, cm.jet)
	x,y=mgrid[1:2,1:10:0.1]
	imshow(y, cmap=light_jet)
	"""
	cdict = cmap._segmentdata
	step_dict = {}
	# Firt get the list of points where the segments start or end
	for key in ('red','green','blue'):
		step_dict[key] = map(lambda x: x[0], cdict[key])
	step_list = sum(step_dict.values(), [])
	step_list = array(list(set(step_list)))
	# Then compute the LUT, and apply the function to the LUT
	reduced_cmap = lambda step : array(cmap(step)[0:3])
	old_LUT = array(map( reduced_cmap, step_list))
	new_LUT = array(map( function, old_LUT))
	# Now try to make a minimal segment definition of the new LUT
	cdict = {}
	for i,key in enumerate(('red','green','blue')):
		this_cdict = {}
		for j,step in enumerate(step_list):
			if step in step_dict[key]:
				this_cdict[step] = new_LUT[j,i]
			elif new_LUT[j,i]!=old_LUT[j,i]:
				this_cdict[step] = new_LUT[j,i]
	colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
	colorvector.sort()
	cdict[key] = colorvector

	return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_discretize(cmap, N):
	"""Return a discrete colormap from the continuous colormap cmap.

	cmap: colormap instance, eg. cm.jet. 
	N: Number of colors.

	Example
	x = resize(arange(100), (5,100))
	djet = cmap_discretize(cm.jet, 5)
	imshow(x, cmap=djet)
	"""

	cdict = cmap._segmentdata.copy()
	# N colors
	colors_i = linspace(0,1.,N)
	# N+1 indices
	indices = linspace(0,1.,N+1)
	for key in ('red','green','blue'):
		# Find the N colors
		D = array(cdict[key])
		I = interpolate.interp1d(D[:,0], D[:,1])
		colors = I(colors_i)
		# Place these colors at the correct indices.
		A = zeros((N+1,3), float)
		A[:,0] = indices
		A[1:,1] = colors
		A[:-1,2] = colors
		# Create a tuple for the dictionary.
		L = []
		for l in A:
			L.append(tuple(l))
		cdict[key] = tuple(L)
	# Return colormap object.
	return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_xmap(function,cmap):
	""" Applies function, on the indices of colormap cmap. Beware, function
	should map the [0, 1] segment to itself, or you are in for surprises.

	See also cmap_xmap.
	"""
	cdict = cmap._segmentdata
	function_to_map = lambda x : (function(x[0]), x[1], x[2])
	for key in ('red','green','blue'):
		cdict[key] = map(function_to_map, cdict[key])
		cdict[key].sort()
		assert (cdict[key][0]<0 or cdict[key][-1]>1), "Resulting indices extend out of the [0, 1] segment."


	return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

