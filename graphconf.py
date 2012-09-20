import traits.api as eta
import traitsui.api as etua
import traitsui.menu as etum
import traitsui.wx.tree_editor as etuwt
import pyface.api as epa
import os
import glob
import math
import pprint
import simplejson as sj



def import_emm_module():
	import mayavi.mlab as emm



################################################################################
#                                                                              #
# Dictionnaries                                                                #
#                                                                              #
################################################################################



color_corr = {
	"b": "blue",
	"g": "green",
	"r": "red",
	"c": "cyan",
	"m": "magenta",
	"y": "yellow",
	"w": "white",
	"k": "black"
}



default_colors = [
	"b",
	"g",
	"r",
	"c",
	"m",
	"y",
	"k"
]



web_colors = [
	'maroon',
	'red',
	'purple',
	'fuchsia',
	'green',
	'lime',
	'olive',
	'yellow',
	'navy',
	'blue',
	'teal',
	'aqua'
]



brown_colors = [
	'cornsilk',
	'blanchedalmond',
	'bisque',
	'navajowhite',
	'wheat',
	'burlywood',
	'tan',
	'rosybrown',
	'sandybrown',
	'goldenrod',
	'darkgoldenrod',
	'peru',
	'chocolate',
	'saddlebrown',
	'sienna',
	'brown',
	'maroon'
]



fancy_colors = [
	"orchid",
	"firebrick",
	"orange",
	"peru",
	"olive",
	"saddlebrown",
	"mediumseagreen",
	"darkgreen",
	"slategrey",
	"skyblue",
	"royalblue"
]



curve_line = {
	"none":"",
	"solid": "-",
	"dash": "--",
	"dot": "--",
	"dash dot": "--",
	"dash dot dot": "--"
}



curve_marker = {
	"none": "",
	"plus": "+",
	"cross": "x",
	"asterisk": "h",
	"circle": ",",
	"square": "s",
	"diamond": "1",
	"triangle up": "2",
	"triangle down": "3",
	"star": "4"
}



################################################################################
#                                                                              #
# Usefull tools.                                                               #
# * colorize                                                                   #
# * get_data                                                                   #
# * rotate_data( |_by_mouse)                                                   #
# * translate_data( | _by_mouse)                                               #
# * get_index_list                                                             #
# * remove_line                                                                #
#                                                                              #
################################################################################



def cf():

	return gcf()



def cAs():

	return getp(cf(), 'axes')



def cA(num = -1):

	return getp(cf(), 'axes')[num]



def ci(image = -1):

	return gci()



def clv():

	to_return = []

	for i in cA().lines:
		if getp(i, 'visible'):
			to_return.append(i)

	return to_return



def cls():

	return clv()



def cv(line_number = -1):

	return clv()[line_number]



def cl(line_number = -1):

	return cv(line_number)



def clh():

	to_return = []

	for i in cA().lines:
		if getp(i, 'visible') == False:
			to_return.append(i)

	return to_return



def ch(line_number = -1):

	return clh()[line_number]



def cla():

	return cA().lines



def ca(line_number = -1):

	return cla()[line_number]



def hide_or_not(i):

	if getp(i, 'visible') == True:
		return ' '

	else:
		return 'H'



def get_info():

	n = 0

	for i in cla():
		print(hide_or_not(i) + " " + str(n)
		+ "\tlen: " + str(len(getp(i, 'xdata')))
		+ "\tmarker: " + getp(i, 'marker')
		+ "\tls: " + getp(i, 'ls')
		+ "\tc: " + getp(i, 'c'))
		n += 1



def ct(text_number = -1):

	return gca().texts[text_number]



def ct(text_number = -1):

	return gca().texts



def cx():
	
	return cA().xaxis



def cy():
	
	return cA().yaxis



def cT():
	
	return ca().title



def cPs():

	return filter(lambda i: type(i) == Polygon, getp(ca(), 'children'))



def cP(num = -1):

	return filter(lambda i: type(i) == Polygon, getp(ca(), 'children'))[num]



def cEs():

	return filter(lambda i: type(i) == Ellipse, getp(ca(), 'children'))



def cE(num = -1):

	return filter(lambda i: type(i) == Ellipse, getp(ca(), 'children'))[num]



def colorize(palette = "fancy", offset = 0, period = None, reverse = False):

	if type(palette) == str:
		if palette == "fancy":
			colors = fancy_colors

		elif palette == "web":
			colors = web_colors

		elif palette == "brown":
			colors = brown_colors

		else:
			colors = default_colors

	elif type(palette) == list:
		colors = palette

	else:
		colors = default_colors

	if period == None:
		period = len(colors)

	n = len(clv())

	if (reverse):
		colors.reverse()

	for i in range(n):
		tmp_col = colors[i%period + offset%(len(colors) - 1)]
		setp(cv(i), 'color', tmp_col, 'mfc', tmp_col, 'mec', tmp_col)

	draw()



def get_data(layer = -1):

	return gca().lines[layer].get_data()



def rotate_data(theta = 0, layer = None):

    if layer == None:
        layer = range(len(gca().lines))

    for j in layer:
        i = gca().lines[j]
        _x1, _y1 = i.get_data()
        i.set_xdata(cos(theta) * _x1 + sin(theta) * _y1)
        i.set_ydata(-sin(theta) * _x1 + cos(theta) * _y1)

    draw()
	


def translate_data(displacement = 0):

	for i in gca().lines:
		_x1, _y1 = i.get_data()
		i.set_xdata(_x1 + displacement[0])
		i.set_ydata(_y1 + displacement[1])

	draw()



def rotate_data_by_mouse():

	_p = ginput(2, show_clicks=False)

	if len(_p) < 2:
		return

	dtheta = - arctan2(_p[1][1], _p[1][0]) + arctan2(_p[0][1], _p[0][0])
	print("dtheta: " + str(dtheta))

	for i in gca().lines:
		_x1, _y1 = i.get_data()
		i.set_xdata(cos(dtheta) * _x1 + sin(dtheta) * _y1)
		i.set_ydata(-sin(dtheta) * _x1 + cos(dtheta) * _y1)

	draw()



def translate_data_by_mouse(layer = None):

	_p = ginput(2, show_clicks=False)

	if len(_p) < 2:
		return

	xdisplacement = _p[1][0] - _p[0][0]
	ydisplacement = _p[1][1] - _p[0][1]
	print("dx: " + str(xdisplacement), "dy :" + str(ydisplacement))

	if layer != None:
		_x1, _y1 = gca().lines[layer].get_data()
		gca().lines[layer].set_xdata(_x1 + xdisplacement)
		gca().lines[layer].set_ydata(_y1 + ydisplacement)	
	
	else:
		for i in gca().lines:
			_x1, _y1 = i.get_data()
			i.set_xdata(_x1 + xdisplacement)
			i.set_ydata(_y1 + ydisplacement)

	draw()



def get_index_list():

	n = 0

	for i in gca().lines:
		print("# " + str(n)+ "\tmarker: " + str(i.get_marker()) + "\tline: " + str(i.get_linestyle()) + "\tcolor: " + str(i.get_color()) + "\tlen: " + str(len(i.get_xdata())))
		n += 1
	
	return len(gca().lines)



def slice_graph(num = None):

	source = cA().lines
	target = figure()

	for i in num:
		x, y = getp(source[i], 'data')
		marker = getp(source[i], 'marker')
		ms = getp(source[i], 'ms')
		ls = getp(source[i], 'ls')
		lw = getp(source[i], 'lw')
		c = getp(source[i], 'c')
		label = getp(source[i], 'label')
		plot(x, y, marker=marker, markersize=ms, linestyle=ls, linewidth=lw, color=c, label=label)



def explode_graph(num = None):

	source = gcf()

	for num in range(len(source.axes[-1].lines)):
		x, y = source.axes[-1].lines[num].get_data()
		marker = source.axes[-1].lines[num].get_marker()
		markersize = source.axes[-1].lines[num].get_markersize()
		linestyle = source.axes[-1].lines[num].get_linestyle()
		linewidth = source.axes[-1].lines[num].get_linewidth()
		color = source.axes[-1].lines[num].get_color()
		label = source.axes[-1].lines[num].get_label()
		figure()
		title("curve: " + str(num) + ", visible: " + str(getp(source.axes[-1].lines[num], 'visible')))
		plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)



def copy_paste_graph(source_fig, drop_fig, xoffset = 0, xmul = 1, yoffset = 0, ymul = 1):

	source = figure(source_fig)
	to_process = clv()
	drop = figure(drop_fig)

	for i in to_process:
		x, y = getp(i, 'data')
		marker = getp(i, 'marker')
		ms = getp(i, 'ms')
		mec = getp(i, 'mec')
		mfc = getp(i, 'mfc')
		ls = getp(i, 'ls')
		lw = getp(i, 'lw')
		c = getp(i, 'c')
		label = getp(i, 'label')
		figure(drop.number)
		plot(xmul*array(x) + xoffset, ymul*array(y) + yoffset)
		setp(cv(), 'marker', marker, 'ms', ms, 'mec', mec, 'mfc', mfc, 'ls', ls, 'lw', lw, 'c', c, 'label', label)



def get_norm_angle():

	a = ginput(show_clicks=False)[0]

	return {
		"norm": sqrt(a[0]**2 + a[1]**2),
		"angle_rad": arctan2(a[1], a[0]),
		"angle_deg": arctan2(a[1], a[0])*180./pi
	}



def get_interval():

	a = ginput(2, show_clicks=False)
	dx = a[1][0] - a[0][0]
	dy = a[1][1] - a[0][1]

	return {
		"interval x": dx,
		"interval y": dy,
		"interval r": sqrt(dx**2 + dy**2)
	}



def remove_first(layer = 0):

	setp(ca(layer), 'visible', False)



def remove_line(layer = None):

	if layer == None:
		layer = [-1]
	
	if type(layer) == int:
		layer = [layer]

	for i in layer:
		setp(ca(i), 'visible', False)

			

def keep_line(layer2keep = -1):
	
	if type(layer2keep) == int:
		layer2keep = [layer2keep]

	for i in range(len(cla())):
		if i in layer2keep:
			continue
		setp(ca(i), 'visible', False)

	

def remove_last(layer = None, axis = -1):

	if layer == None:
		setp(ca(-1), 'visible', False)
	else:
		setp(ca(layer))
	draw()



def hide_line(layer = -1):

	if type(layer) == int:
		layer = [layer]

	for i in layer:
		setp(ca(i), 'visible', False)



def show_line(layer = -1):

	if type(layer) == int:
		layer = [layer]

	for i in layer:
		setp(ca(i), 'visible', True)



def show_all_lines():

    for i in cla():
        setp(i, 'visible', True)



def hide_all_lines():

    for i in cla():
        setp(i, 'visible', False)



def remove_last_patch():

	if len(cA().patches) > 0:
		cA().patches.pop()
		draw()



def fastgen(a, b, cond=True):

	return (i for i in range(a, b) if cond)



def line_matrix(xstep = 1, ystep = 1, bottom = None, top = None):

	xmin = getp(cv(0), 'xdata').min()
	xmax = getp(cv(0), 'xdata').max()
	ymin = getp(cv(0), 'ydata').min()
	ymax = getp(cv(0), 'ydata').max()

	xx, yy = meshgrid(arange(xmin, xmax, xstep), arange(ymin, ymax, ystep))

	cc = 0.0*xx

	print("grid xx: " + str(len(xx)) + "x" + str(len(xx[0])))
	print("grid yy: " + str(len(yy)) + "x" + str(len(yy[0])))
	print("grid cc: " + str(len(cc)) + "x" + str(len(cc[0])))

	intx = xrange(0, int(round((xmax - xmin)/xstep)))
	inty = xrange(0, int(round((ymax - ymin)/ystep)))

	for j in intx:
		for k in inty:
			cc[k][j] = 0

	if len(clv()) == 1:
		offset = 1
	else:
		offset = 0

	for n in xrange(len(clv())):
		if bottom != None:
			if n < bottom:
				n = bottom
		if top != None:
			if n > top:
				n = top
		data = getp(cv(n), 'xydata')
		for i in xrange(len(data)):
				j = (data[i][1] - ymin)/ystep
				k = (data[i][0] - xmin)/xstep
				# print("j, k: " + str(j) + ", " + str(k))
				try:
					cc[round(j)][round(k)] = n + offset
				except:
					None #print("oups, stopped @ " + str(i) + "/" + str(len(data)))

	return xx, yy, cc



def sort_filename(list2sort):

	N = len(list2sort)

	good = []
	bad = []

	for n in range(N):
		try:
			if list2sort[n][-3:-1] == 'lo':
				good.append(list2sort[n])
			else:
				bad.append(list2sort[n])
		except:
			print("stuck", n)
        
	return {'good': good, 'bad': bad}



def process_plot_to_png():

	l = sort_filename(os.listdir("."))

	for i in l["good"]:
		openGraph(i)
		if len(clv()) > 0: 
			savefig(i + ".png")
		close()



def process_plot_to_pdf():

	l = sort_filename(os.listdir("."))

	for i in l["good"]:
		openGraph(i)
		if len(clv()) > 0: 
			savefig(i + ".pdf")
		close()



################################################################################
#                                                                              #
# Different class of usefull graph.                                            #
# * graph(...)                                                                 #
# * (im | p)graph(...)                                                         #
#                                                                              #
################################################################################



import pylab



def tight_axis():

	xmin = Inf
	xmax = -Inf
	ymin = Inf
	ymax = -Inf
	
	for i in clv():
		x, y = getp(i, 'data')
		if x.min() < xmin and x.min() != -Inf:
			xmin = x.min()
		if x.max() > xmax and x.max() != Inf:
			xmin = x.max()
		if y.min() < ymin and y.min() != -Inf:
			ymin = y.min()
		if y.max() < ymax and y.max() != Inf:
			ymax = y.max()
	
	if xmin != Inf:
		axis(xmin=xmin)
	
	if xmax != -Inf:
		axis(xmax=xmax)
	
	if ymin != Inf:
		axis(ymin=ymin)
	
	if ymax != -Inf:
		axis(ymax=ymax)
	
	print (xmin, xmax, ymin, ymax)



def show_colormap():
	
	a = outer(arange(0, 1, 0.01), ones(10))
	figure(figsize = (10, 5))
	subplots_adjust(top = 0.8, bottom = 0.05, left = 0.01, right = 0.99)
	maps=[m for m in cm.datad if not m.endswith("_r")]
	maps.sort()
	l = len(maps) + 1

	for i, m in enumerate(maps):
		subplot(1, l, i + 1)
		axis("off")
		imshow(a, aspect='auto', cmap=get_cmap(m), origin="lower")
		title(m, rotation=45, fontsize=10)



def graph(*args):
	
	_f = figure()
	plot(*args)
	return _f



class imgraph:

	def __init__(self, *args):

		figure()
		imshow(*args)



class pgraph:

	def __init__(self, *args, **kwargs):

		figure()
		pcolormesh(*args, **kwargs)



def contourgraph(*args, **kwargs):

	f = figure()
	contour(*args, **kwargs)
	return f



def contourfgraph(*args, **kwargs):

	f = figure()
	contourf(*args, **kwargs)
	return f



class polargraph:

	def __init__(self, *args, **kwargs):

		figure()
		polar()
		pcolormesh(*args, **kwargs)



def savepdfpng(fname, dpi = 300):

	print "save png"
	setp(gci(), 'visible', True)
	axis('off')
	savefig(fname + '.png', dpi = dpi)
	print "save pdf"
	setp(gci(), 'visible', False)
	axis('on')
	savefig(fname + '.pdf', dpi = dpi)
	setp(gci(), 'visible', True)
	axis('on')
	print "done"



class getLasso:

	def __init__(self):

		self.xs = []
		self.ys = []
		self.cid = cf().canvas.mpl_connect('button_press_event', self._click)
		self.kid = cf().canvas.mpl_connect('key_press_event', self._press)

	def _click(self, event):

		if event.inaxes != cA():
			return

		self.xs.append(event.xdata)
		self.ys.append(event.ydata)

		if len(self.xs) == 1:
			cA().add_patch(Polygon(zip(self.xs, self.ys), True, alpha = 0.4))
			draw()

		if len(self.xs) > 1:
			setp(cA().patches[-1], 'xy', zip(self.xs, self.ys))

	def _press(self, event):

		remove_last_patch()
		pprint.pprint(zip(self.xs, self.ys))
		cf().canvas.mpl_disconnect(self.cid)
		cf().canvas.mpl_disconnect(self.kid)
		draw()
	
	def stop(self):

		cf().canvas.mpl_disconnect(self.cid)
		cf().canvas.mpl_disconnect(self.kid)



class getSquare:

	def __init__(self):

		self.x = []
		self.y = []
		self.a = None
		self.b = None
		self.c = None
		self.d = None
		self.first_click = True
		self.cid = cf().canvas.mpl_connect('button_press_event', self._click)

	def _click(self, event):

		if event.inaxes != cA():
			return

		if self.first_click:
			self.a = event.xdata
			self.b = event.ydata
			self.x.append(self.a)
			self.y.append(self.b)
			self.first_click = False
			cA().add_patch(Polygon(zip(self.x, self.y), True, alpha = 0.4))
			draw()
		else:
			self.c = event.xdata
			self.d = event.ydata

			self.x.append(self.c)
			self.y.append(self.b)
			self.x.append(self.c)
			self.y.append(self.d)
			self.x.append(self.a)
			self.y.append(self.d)
			self.x.append(self.a)
			self.y.append(self.b)
			
			setp(cA().patches[-1], 'xy', zip(self.x, self.y))
			setp(cA().patches[-1], 'xy', zip(self.x, self.y))
			setp(cA().patches[-1], 'xy', zip(self.x, self.y))
			setp(cA().patches[-1], 'xy', zip(self.x, self.y))
			# remove_last_patch()
			pprint.pprint(zip(self.x, self.y))
			cf().canvas.mpl_disconnect(self.cid)
			draw()
	
	def stop(self):

		cf().canvas.mpl_disconnect(self.cid)



class LassoFilter:

	def __init__(self):

		self.inside = True
		self.xs = []
		self.ys = []
		self.cid = cf().canvas.mpl_connect('button_press_event', self._click)
		self.kid = cf().canvas.mpl_connect('key_press_event', self._press)

	def _click(self, event):

		if event.inaxes != cA():
			return

		self.xs.append(event.xdata)
		self.ys.append(event.ydata)

		if len(self.xs) == 1:
			cA().add_patch(Polygon(zip(self.xs, self.ys), True, alpha = 0.4))
			draw()

		if len(self.xs) > 1:
			setp(cA().patches[-1], 'xy', zip(self.xs, self.ys))

	def _press(self, event):

		if event.key == 'x':
			print("take outside into account")
			self.inside = False
		else:
			print("take inside into account")
			self.inside = True

		remove_last_patch()
		lasso_lines(zip(self.xs, self.ys), self.inside)
		cf().canvas.mpl_disconnect(self.cid)
		cf().canvas.mpl_disconnect(self.kid)
		draw()
	
	def stop(self):

		cf().canvas.mpl_disconnect(self.cid)
		cf().canvas.mpl_disconnect(self.kid)



def lasso_filter(lasso, points, inside = True):
	
	verts = array(lasso)
	line_inside = []
	line_outside = []

	for i in range(len(points)):
		if matplotlib.nxutils.pnpoly(points[i][0], points[i][1], verts) == True:
			line_inside.append(points[i].tolist())
		else:
			line_outside.append(points[i].tolist())

	if inside:
		return line_inside
	else:
		return line_outside



def lasso_lines(lasso, inside=True):
	
	lines_to_remove = []

	for i in range(len(cla())):
		new_line = lasso_filter(lasso, getp(ca(i), 'xydata'), inside)
		x, y = hsplit(array(new_line), 2)

		if len(x) == 0:
			lines_to_remove.insert(0, i)

		setp(ca(i), 'data', [x, y])

	for i in lines_to_remove:
		cA().lines.pop(i)



################################################################################
#                                                                              #
# Importating and saving graph and data.                                       #
# * (save | open)JSONgraph(...)                                                #
# * openNanoQt(graph | polar_graph)(...)                                       #
#                                                                              #
################################################################################


key_words = {
	"figure" : ['alpha', 'animated','edgecolor','facecolor','figheight','figwidth','frameon','label','visible','zorder'],
	"axes" : ['adjustable', 'alpha', 'anchor', 'animated', 'aspect', 'autoscale_on', 'axis_bgcolor', 'axisbelow', 'frame_on', 'label', 'legend', 'title', 'visible', 'xbound', 'xlabel', 'xscale', 'ybound', 'ylabel', 'yscale', 'zorder'],
	"lines" : ['alpha', 'animated', 'antialiased', 'color', 'dash_capstyle', 'dash_joinstyle', 'drawstyle', 'label', 'linestyle', 'linewidth', 'marker', 'markeredgecolor', 'markeredgewidth', 'markerfacecolor', 'markersize', 'solid_capstyle', 'solid_joinstyle', 'visible', 'zorder'],
	"texts" : ['alpha', 'animated', 'color', 'family', 'horizontalalignment', 'label', 'position', 'rotation', 'size', 'stretch', 'style', 'text', 'variant', 'verticalalignment', 'visible', 'weight', 'zorder']
}


def saveJSONgraph(file_name = None):

	if file_name == None:
		print("Filename needed")
		return

	# Here is the bundle we plan to save
	bundle = {}

	# First, save figure parameters coming from gcf()
	bundle['figure'] = {} # place to save figure's properties

	for i in key_words["figure"] :
		bundle["figure"][i] = getp(gcf(), i)

	# Now, save axes.
	# Axes composed with properties with lines (data) and texts (lables).
	# Need a loop to scan all of the stuffs in all of the axes.
	bundle['axes'] = [] # list of axes

	for current_axis in gcf().axes:
		# First of all, pack the general properties of the axes.
		bundle['axes'].append({})
		
		for i in key_words["axes"] :
			bundle['axes'][-1][i] = getp(current_axis, i)
		bundle['axes'][-1]['position'] = [current_axis._position._get_x0(), current_axis._position._get_y0(), current_axis._position._get_width(), current_axis._position._get_height()]

		# Now, pack the line.
		bundle['axes'][-1]['lines'] = []

		for current_line in current_axis.lines:
			bundle['axes'][-1]['lines'].append({})

			for j in key_words["lines"] :
				bundle['axes'][-1]['lines'][-1][j] = getp(current_line, j)
			bundle['axes'][-1]['lines'][-1]['xdata'] = getp(current_line, 'xdata').tolist()
			bundle['axes'][-1]['lines'][-1]['ydata'] = getp(current_line, 'ydata').tolist()

		# Then, pack the texts
		bundle['axes'][-1]['texts'] = []

		for current_text in current_axis.texts:
			bundle['axes'][-1]['texts'].append({})
			
			for j in key_words["texts"] :
				bundle['axes'][-1]['texts'][-1][j] = getp(current_text, j)

	target_file = file(file_name, "w")
	sj.dump(bundle, target_file)
	target_file.close()



def openJSONgraph(bundle = None):

	try:
		setp(gcf(), 'alpha', bundle['figure']['alpha'], 
			    'animated', bundle['figure']['animated'], 
			    'edgecolor', bundle['figure']['edgecolor'], 
			    'facecolor', bundle['figure']['facecolor'], 
			    'frameon', bundle['figure']['frameon'], 
			    'label', bundle['figure']['label'], 
			    'visible', bundle['figure']['visible'], 
			    'zorder', bundle['figure']['zorder']) 
	except:
		print("Warning: fail to set figure properties")

	# Now, set axes.
	# Axes composed with properties with lines (data) and texts (lables).
	# Need a loop to scan all of the stuffs in all of the axes.

	for current_axis in bundle['axes']:
		for current_line in current_axis['lines']:
			try:
				plot(current_line['xdata'], current_line['ydata'])	
			except:
				print("Warning: fail to plot data")

			try:
				setp(gca().lines[-1], 'alpha'  , current_line['alpha']          ,
						      'animated'       , current_line['animated']       ,
						      'antialiased'    , current_line['antialiased']    ,
						      'color'          , current_line['color']          ,
						      'dash_capstyle'  , current_line['dash_capstyle']  ,
						      'dash_joinstyle' , current_line['dash_joinstyle'] ,
						      'drawstyle'      , current_line['drawstyle']      ,
						      'label'          , current_line['label']          ,
						      'linestyle'      , current_line['linestyle']      ,
						      'linewidth'      , current_line['linewidth']      ,
						      'marker'         , current_line['marker']         ,
						      'markeredgecolor', current_line['markeredgecolor'],
						      'markeredgewidth', current_line['markeredgewidth'],
						      'markerfacecolor', current_line['markerfacecolor'],
						      'markersize'     , current_line['markersize']     ,
						      'solid_capstyle' , current_line['solid_capstyle'] ,
						      'solid_joinstyle', current_line['solid_joinstyle'],
						      'visible'        , current_line['visible']        ,
						      'zorder'         , current_line['zorder']         )
			except:
				print("Warning: fail to set line properties")

		for current_text in current_axis['texts']:		
			try:	
				text(current_text['position'][0], current_text['position'][1], current_text['text'])
			except:
				print("Warning: fail to set text")

			try:
				setp(gca().texts[-1], 'alpha'              , current_text['alpha']              , 

						      'animated'           , current_text['animated']           ,
						      'color'              , current_text['color']              ,
						      'family'             , current_text['family']             ,
						      'horizontalalignment', current_text['horizontalalignment'],
						      'label'              , current_text['label']              ,
						      'rotation'           , current_text['rotation']           ,
						      'size'               , current_text['size']               ,
						      'stretch'            , current_text['stretch']            ,
						      'style'              , current_text['style']              ,
						      'variant'            , current_text['variant']            ,
						      'verticalalignment'  , current_text['verticalalignment']  ,
						      'visible'            , current_text['visible']            ,
						      'weight'             , current_text['weight']             ,
						      'zorder'             , current_text['zorder']             )
			except:
				print("Warning: fail to set text properties")

		try:
			setp(gca(), 'adjustable'  , current_axis['adjustable']  ,
				    'alpha'       , current_axis['alpha']       ,
				    'anchor'      , current_axis['anchor']      ,
				    'animated'    , current_axis['animated']    ,
				    'aspect'      , current_axis['aspect']      ,
				    'autoscale_on', current_axis['autoscale_on'],
				    'axis_bgcolor', current_axis['axis_bgcolor'],
				    'axisbelow'   , current_axis['axisbelow']   ,
				    'frame_on'    , current_axis['frame_on']    ,
				    'label'       , current_axis['label']       ,
				    'position'    , current_axis['position']    ,
				    'title'       , current_axis['title']       ,
				    'visible'     , current_axis['visible']     ,
				    'xlabel'      , current_axis['xlabel']      ,
				    'xscale'      , current_axis['xscale']      ,
				    'ylabel'      , current_axis['ylabel']      ,
				    'yscale'      , current_axis['yscale']      ,
				    'zorder'      , current_axis['zorder']      )
		except:
			print("Warning: fail to set axes properties")

		if type(current_axis['xbound']) == list and type(current_axis['ybound']) == list:
			try:
				setp(gca(), 'xbound'      , current_axis['xbound']      ,
					    'ybound'      , current_axis['ybound']      )
			except:
				print("Warning: Get stuck with (x|y)bounds")



def openNanoQtgraph(bundle, reversed = False):

	debug_string = ""

	try :
		print('Set the title')
		title(bundle['title'])
	except:
		print('Warning: Fail to set the title: ' + bundle['title'])

	try :
		print('Set the x label')
		xlabel(bundle['x_label'])
	except:
		print('Warning: Fail to set the x label: ' + bundle['x_label'])

	try :
		print('Set the y label')
		ylabel(bundle['x_label'])
	except:
		print('Warning: Fail to set the y label: ' + bundle['y_label'])

	try: 
		print('Set the xscale')

		if bundle['logscale_x']:
			setp(cA(), 'xscale', 'log')
		else:
			setp(cA(), 'xscale', 'linear')
	except: 
		print('Warning: fail to set the xscale')

	try: 
		print('Set the yscale')
		if bundle['logscale_y']:
			setp(cA(), 'yscale', 'log')
		else:
			setp(cA(), 'yscale', 'linear')
	except: 
		print('Warning: Fail to set the yscale')

	try:
		print("Plot lines...")

		for n in range(len(bundle["curves"])):
			if bundle["curves"][n] != None:
				
				curvex, curvey = hsplit(array(bundle["curves"][n]["data"]), 2)
				plot(curvex, curvey)

				try:
					print("Set color")
			
					if str(bundle["curves"][n]["options"]["pen_color"]).startswith("0x"):
						formatted_color = '#' + str(bundle["curves"][n]["options"]["pen_color"])[2:]
					else:
						formatted_color = str(bundle["curves"][n]["options"]["pen_color"])
						setp(ca(), 'c', formatted_color)
				except:
					formatted_color = getp(ca(), 'c')
					print("Warning: Fail to set color")

				setp(ca(), 'mec', formatted_color)
				setp(ca(), 'mfc', formatted_color)

				try:
					print("Set line style")
					pen_style = bundle["curves"][n]["options"]["pen_style"]
					line_style = curve_line[pen_style]
					setp(ca(), 'ls', line_style)
				except:
					setp(ca(), 'ls', '')
					print("Warning: Fail to set line style")

				try:
					print("Set marker")
					setp(ca(), 'marker', curve_marker[bundle["curves"][n]["options"]["symbol"]])
				except:
					setp(ca(), 'marker', 's')
					print("Warning: Fail to set marker")

				try:
					print("Set markersize")
					setp(ca(), 'markersize', bundle["curves"][n]["options"]["symbol_size"])

				except:
					setp(ca(), 'markersize', 1)
					print("Warning: Fail to set markersize")

				try:
					print("Set axis")
					axis((bundle["x_min"], bundle["x_max"], bundle["y_min"], bundle["y_max"]))

				except:
					print("Warning: Fail to set axis")
					
	except:
		print("Warning: Fail to plot lines...")

	try:
		print("Draw array...")
		a = array(bundle["array"]["data"]) 
		x, y = meshgrid(linspace(bundle["x_min"], bundle["x_max"], len(a)), linspace(bundle["y_min"], bundle["y_max"], len(a[0])))

		if reversed:
			pcolormesh(x, y, a.max() - a)

		else:
			pcolormesh(x, y, a)

	except:
		print("Warning: Fail to draw array...")

	return



def openNanoQtgraph_polar(bundle, reversed = False, deg_rad = 1):

	print("openNanoQtgraph_polar")
    
	polar()

	try:
		print("Draw array...")
		a = array(bundle["array"]["data"]) 
		x, y = meshgrid(linspace(bundle["x_min"], bundle["x_max"], len(a)), linspace(bundle["y_min"], bundle["y_max"], len(a[0])))

		if reversed:
			pcolormesh(deg_rad*y, x, a.max() - a)
		else:
			pcolormesh(deg_rad*y, x, a)

	except:
		print("Warning: Fail to draw array...")



def openGraph(fname, reversed = False, tight=0):

	if fname == "last":
		_l = os.listdir(".")
		_l.sort()
		fname = [_l.pop()]
		print(fname)
	else:
		_l = glob.glob(fname)
		print(_l)

		if len(_l) < 1:
			print("don't find occurence")
			return

		if len(_l) > 1:
			print("many occurences")
			print("would you like to open everything?")
			print("(y/n)")

			if raw_input() != "y":
			    return

		fname = _l

	for i in fname:
		print(i)

		try:
			print("Open the file")
			f = file(i, "r")
		except:
			print("Warning: Fail to open the file")
			return

		try:
			print("Deseriliaze the JSON string")
			data = sj.load(f)
		except:
			print("Warning: Fail to deseriliaze the JSON string")
			return

		_f = figure ()

		try:
			print("Open the graph in Python format")
			openJSONgraph(data)

			if tight == 1:
				axis("tight")
		except:
			print("Warning: Fail to open the graph in Python format")

		try:
			print("Open the graph in NanoQt format!")
			openNanoQtgraph(data, reversed = reversed)
			if tight == 1:
				axis("tight")
		except:
			print("Warning: Fail to open the graph in NanoQt format!")



def openPolar(fname, reversed = False, angle = 2 * pi):

	if fname == "last":
		_l = os.listdir(".")
		_l.sort()
		fname = _l.pop()
		print(fname)
	else:
		_l = glob.glob(fname)

		if (len(_l) > 1):
			print("too much occurence")
			print(_l)
			return

		if len(_l) < 1:
			print("don't find occurence")
			return

		fname = _l.pop()
	print(fname)

	try:
		print("Open the file")
		f = file(fname, "r")
	except:
		print("Warning: Fail to open the file")
		return

	try:
		print("Deseriliaze the JSON string")
		data = sj.load(f)
	except:
		print("Warning: Fail to open Polar plot")
		return

	figure()

	try:
		print("Processing")
		# openNanoQtgraph_polar(data, reversed = reversed, deg_rad = eval(angle))
		openNanoQtgraph_polar(data, reversed, angle)
	except:
		print("Warning: ouch @ processing!")



class UiImport(eta.HasTraits):

    data = eta.Dict()
    fname = eta.File()
    angle = eta.Str("2*pi")
    reversed = eta.Bool()
    open_graph = eta.Button("Open graph")
    open_polar = eta.Button("Open polar")
 
    def __init__(self, fname = ""):

        self.fname = os.getcwdu()
        self.directory = os.getcwdu()
        self.current_directory = self.directory

    def _open_graph_fired(self):

        openGraph(self.fname, self.reversed)

    def _open_polar_fired(self):

        try:
            self.data = sj.load(file(self.fname, "r"))

            try:
                openNanoQtgraph_polar(self.data, reversed = self.reversed, deg_rad = eval(self.angle))

            except:
                print("Warning: ouch @ processing!")

        except:
            print("Warning: ouch @ opening!")

    open_graph = etum.Action(name = 'open graph', action = '_open_graph_fired')

    open_polar = etum.Action(name = 'open polar', action = '_open_polar_fired')

    view = etua.View(
        etua.Item('angle'),
        etua.Item('reversed'),
        etua.Item('fname', editor=etua.FileEditor(auto_set = True), style = "custom"),
        resizable = True,
        scrollable = True,
		title= "UiImport",
		toolbar = etum.ToolBar(open_graph, open_polar),
        height = 720,
        width = 800
    )



def uiimport(fname = ""):

    wxuiimport = UiImport(fname)
    wxuiimport.configure_traits()



class UiBrowser(eta.HasTraits):

	directory = eta.Directory()
	current_directory = eta.Str("")
	cd = eta.Button("Directory")
	plot_to_png = eta.Button(".plot to .png")

	def __init__(self):

		self.directory = os.getcwdu()
		self.current_directory = self.directory

	def _cd_fired(self):

		_ip.magic("cd " + self.directory)
		self.current_directory = os.getcwdu() # _ip.magic("os.getcwdu() ")

	def _plot_to_png_fired(self):

		self.current_directory = os.getcwdu() # _ip.magic("pwd ")
		print("Convert .plot to .png in " + os.getcwdu()) # _ip.magic("pwd"))
		process_plot_to_png()
		print("Convert finished")

	change_directory = etum.Action(name = 'Directory', action = '_cd_fired')
	convert_plot_png = etum.Action(name = '.plot to .png', action = '_plot_to_png_fired')

	view = etua.View(
		etua.Item('current_directory'),
		etua.Item('directory', editor=etua.DirectoryEditor(), style = "custom"),
		resizable = True,
		scrollable = True,
		title= "Browser",
		toolbar = etum.ToolBar(change_directory, convert_plot_png),
		height = 640,
		width = 800
	)



def browser():

    wxbrowser = UiBrowser()
    wxbrowser.configure_traits()



class Inbox(eta.HasTraits):

    inbox = eta.List
   
    def __init__(self, inbox):

		self.inbox = [str(msg) for msg in inbox]

    view = etua.View(
        etua.Item('inbox', style='custom'),
        resizable = True,
        scrollable = True,
		title= "Inbox",
        height = 640,
        width = 800
    )



def inbox(inb):

    wxinbox = Inbox(inb)
    wxinbox.configure_traits()
