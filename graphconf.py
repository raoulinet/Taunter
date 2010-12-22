import enthought.traits.api as eta
import enthought.traits.ui.api as etua
import simplejson as sj


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
	"none":"None",
	"solid": "-",
	"dash": "--",
	"dot": "--",
	"dash dot": "--",
	"dash dot dot": "--"
	}


curve_marker = {
	"none": "None",
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

def strlist(a):
	"""
	convert a tuple, array or any list like to a string list
	keep comma in the list representation
	"""

	return str(list(a))


def hex_diese_to_tuple(a):
	"""
	Convert a hexa color #8B0000 to RGB tuple
	(255, 255, 255)
	"""

	a.replace("#", "")
	return (eval("0x" + a[0:2]), eval("0x" + a[2:4]), eval("0x" + a[4:6]))


def hex_x_to_tuple(a):
	"""
	Convert a hexa color 0x8B0000 to RGB tuple
	(255, 255, 255)
	"""

	a.replace("0x", "")
	return (eval("0x" + a[0:2]), eval("0x" + a[2:4]), eval("0x" + a[4:6]))


def hex_diese_to_normalized_tuple(a):
	"""
	Convert a hexa color #8B0000 to normalized RGB tuple
	(1.0, 1.0, 1.0)
	"""

	a.replace("#", "")
	return (eval("0x" + a[4:6])/255., eval("0x" + a[2:4])/255., eval("0x" + a[0:2])/255.)


def hex_x_to_normalized_tuple(a):
	"""
	Convert a hexa color #8B0000 to normalized RGB tuple
	(1.0, 1.0, 1.0)
	"""

	a.replace("0x", "")
	return (eval("0x" + a[4:6])/255., eval("0x" + a[2:4])/255., eval("0x" + a[0:2])/255.)


def colorize(palette = "fancy", offset = 0):
	"""
	to colorize lines depdening of layers
	"""
	if palette == "fancy":
		colors = fancy_colors
	else:
		if palette == "web":
			colors = web_colors
		else:
			if palette == "brown":
				colors = brown_colors
			else:
				colors = default_colors

	n = len(gca().lines)
	for i in range(n):
		tmp_col = colors[(i + offset)%len(colors)]
		gca().lines[i].set_color(tmp_col)
		gca().lines[i].set_markerfacecolor(tmp_col)
		gca().lines[i].set_markeredgecolor(tmp_col)
	draw()


def get_data(layer = -1):
	return gca().lines[layer].get_data()


def rotate_data(theta = 0):
	for i in gca().lines:
		_x1, _y1 = i.get_data()
		i.set_xdata(cos(theta) * _x1 + sin(theta) * _y1)
		i.set_ydata(-sin(theta) * _x1 + cos(theta) * _y1)
	draw()
	


def move_data(displacement = 0):
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


def move_data_by_mouse(layer = None):
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
	print("layer count: " + str(len(gca().lines)))
	n = 0
	for i in gca().lines:
		print("# " + str(n)+ ", marker: " + str(i.get_marker()) + ", line: " + str(i.get_linestyle()) + ", color: " + str(i.get_color()) + ", len: " + str(len(i.get_xdata())))
		n += 1


def set_label(layer, label):
	gca().lines[layer].set_label(label)
	draw()


def slice_graph(begin = 0, end = None, step = 1):
	source = gcf()
	drop = figure()
	if end == None:
		end = len(source.axes[-1].lines)
	for i in range(begin, end, step):
		x, y = source.axes[-1].lines[i].get_data()
		marker = source.axes[-1].lines[i].get_marker()
		markersize = source.axes[-1].lines[i].get_markersize()
		linestyle = source.axes[-1].lines[i].get_linestyle()
		linewidth = source.axes[-1].lines[i].get_linewidth()
		color = source.axes[-1].lines[i].get_color()
		label = source.axes[-1].lines[i].get_label()
		figure(drop.number)
		plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)	


def copy_paste_graph(source_fig, drop_fig):
	source = figure(source_fig)
	drop = figure(drop_fig)
	for i in range(len(source.axes[-1].lines)):
		x, y = source.axes[-1].lines[i].get_data()
		marker = source.axes[-1].lines[i].get_marker()
		markersize = source.axes[-1].lines[i].get_markersize()
		linestyle = source.axes[-1].lines[i].get_linestyle()
		linewidth = source.axes[-1].lines[i].get_linewidth()
		color = source.axes[-1].lines[i].get_color()
		label = source.axes[-1].lines[i].get_label()
		figure(drop.number)
		plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)	


def get_norm_angle():
    a = ginput()[0]
    return {
        "norm": sqrt(a[0]**2 + a[1]**2),
        "angle_rad": arctan2(a[1], a[0]),
        "angle_deg": arctan2(a[1], a[0])*180./pi
    }


def remove_line(layer = -1):
	gca().lines.pop(layer)
	draw()


def set_line(layer = None, line = None, width = None, color = None):
    if layer == None:
        for i in gca().lines:
            if line:
                i.set_linestyle(line)
            if width:
                i.set_linewidth(width)
    else:
        if line:
            gca().lines[layer].set_linestyle(line)
        if width:
            gca().lines[layer].set_linewidth(width)
        if color:
            gca().lines[layer].set_color(color)
    draw()


def set_marker(layer = None, marker = None, size = None, color = None):
    if layer == None:
        for i in gca().lines:
            if marker:
                i.set_marker(marker)
            if size:
                i.set_markersize(size)
    else:
        if maker:
            gca().lines[layer].set_marker(marker)
        if size:
            gca().lines[layer].set_markersize(size)
        if color:
            gca().lines[layer].set_markerfacecolor(color)
            gca().lines[layer].set_markeredgecolor(color)
    draw()



def line_matrix(bottom = None, top = None, xstep = 1, ystep = 1):

    xmin = gca().lines[0].get_xdata().min()
    xmax = gca().lines[0].get_xdata().max()
    ymin = gca().lines[0].get_ydata().min()
    ymax = gca().lines[0].get_ydata().max()


    xx, yy = meshgrid(arange(xmin, xmax, xstep), arange(ymin, ymax, ystep))

    cc = 0*xx

    print("grid xx: " + str(len(xx)) + "x" + str(len(xx[0])))
    print("grid yy: " + str(len(yy)) + "x" + str(len(yy[0])))
    print("grid cc: " + str(len(cc)) + "x" + str(len(cc[0])))

    intx = range((xmax - xmin)/xstep)
    inty = range((ymax - ymin)/ystep)

    # print("inty, intx: " + str(len(inty)) + ", " + str(len(intx)))

    for n in range(len(gca().lines)):
        if bottom != None:
            if n < bottom:
                n = bottom
        if top != None:
            if n > top:
                n = top
        data = gca().lines[n].get_xydata()
        for i in range(len(data)):
                j = (data[i][1] - ymin)/ystep
                k = (data[i][0] - xmin)/xstep
                # print("j, k: " + str(j) + ", " + str(k))
                try:
                    cc[round(j)][round(k)] = n
                except:
                    print("oups, stopped @ " + str(i) + "/" + str(len(data)))

    return xx, yy, cc



################################################################################

class Code(eta.HasTraits):
	code = eta.Code
	
	def __init__(self, code = "two words"):
		self.code = code

	def _add_code(self, code):
		_code = self.code
		self.code = _self + code

	traits_view = etua.View(
		etua.Item("code", style="custom"),
		resizable = True,
		scrollable = True,
		width = 640,
		height = 480
	)



################################################################################

rcParams = {
	'agg.path.chunksize': 0,
	'axes.axisbelow': False,
	'axes.edgecolor': 'k',
	'axes.facecolor': 'w',
	'axes.formatter.limits': [-7, 7],
	'axes.grid': False,
	'axes.hold': True,
	'axes.labelcolor': 'k',
	'axes.labelsize': 'x-large', # 
	'axes.linewidth': 1.0,
	'axes.titlesize': 'x-large',
	'axes.unicode_minus': True,
	'backend': 'Agg',
	'backend_fallback': True,
	'cairo.format': 'png',
	'contour.negative_linestyle': 'dashed',
	'datapath': '/Library/Frameworks/Python.framework/Versions/4.3.0/lib/python2.5/site-packages/matplotlib-0.98.5.2n2-py2.5-macosx-10.3-fat.egg/matplotlib/mpl-data',
	'docstring.hardcopy': False,
	'figure.autolayout': False,
	'figure.dpi': 80,
	'figure.edgecolor': 'w',
	'figure.facecolor': '0.75',
	'figure.figsize': [8.0, 6.0],
	'figure.subplot.bottom': 0.10000000000000001,
	'figure.subplot.hspace': 0.20000000000000001,
	'figure.subplot.left': 0.125,
	'figure.subplot.right': 0.90000000000000002,
	'figure.subplot.top': 0.90000000000000002,
	'figure.subplot.wspace': 0.20000000000000001,
	'font.cursive': ['Apple Chancery',
		  'Textile',
		  'Zapf Chancery',
		  'Sand',
		  'cursive'],
	'font.family': 'serif', # 'sans-serif'
	'font.fantasy': ['Comic Sans MS',
		  'Chicago',
		  'Charcoal',
		  'ImpactWestern',
		  'fantasy'],
	'font.monospace': ['Bitstream Vera Sans Mono',
		    'DejaVu Sans Mono',
		    'Andale Mono',
		    'Nimbus Mono L',
		    'Courier New',
		    'Courier',
		    'Fixed',
		    'Terminal',
		    'monospace'],
	'font.sans-serif': ['Bitstream Vera Sans',
		     'DejaVu Sans',
		     'Lucida Grande',
		     'Verdana',
		     'Geneva',
		     'Lucid',
		     'Arial',
		     'Helvetica',
		     'Avant Garde',
		     'sans-serif'],
	'font.serif': ['Bitstream Vera Serif',
		'DejaVu Serif',
		'New Century Schoolbook',
		'Century Schoolbook L',
		'Utopia',
		'ITC Bookman',
		'Bookman',
		'Nimbus Roman No9 L',
		'Times New Roman',
		'Times',
		'Palatino',
		'Charter',
		'serif'],
	'font.size': 20.0, # 12.0
	'font.stretch': 'normal',
	'font.style': 'normal',
	'font.variant': 'normal',
	'font.weight': 'normal',
	'grid.color': 'k',
	'grid.linestyle': ':',
	'grid.linewidth': 1, # 0.5
	'image.aspect': 'equal',
	'image.cmap': 'spectral', # 'jet'
	'image.interpolation': 'nearest', # 'bilinear'
	'image.lut': 256,
	'image.origin': 'lower', # 'upper'
	'image.resample': False,
	'interactive': False,
	'legend.axespad': 0.5,
	'legend.borderaxespad': 0.5,
	'legend.borderpad': 0.40000000000000002,
	'legend.columnspacing': 2.0,
	'legend.fancybox': True, # False
	'legend.fontsize': 'x-large', # large
	'legend.handlelen': 0.050000000000000003,
	'legend.handlelength': 2.0,
	'legend.handletextpad': 0.80000000000000004,
	'legend.handletextsep': 0.02,
	'legend.isaxes': True,
	'legend.labelsep': 0.01,
	'legend.labelspacing': 0.5,
	'legend.loc': 'upper right',
	'legend.markerscale': 1.0,
	'legend.numpoints': 2,
	'legend.pad': 0,
	'legend.shadow': False,
	'lines.antialiased': True,
	'lines.color': 'b',
	'lines.dash_capstyle': 'butt',
	'lines.dash_joinstyle': 'miter',
	'lines.linestyle': '-',
	'lines.linewidth': 1.0,
	'lines.marker': 'None',
	'lines.markeredgewidth': 0.5,
	'lines.markersize': 6,
	'lines.solid_capstyle': 'projecting',
	'lines.solid_joinstyle': 'miter',
	'maskedarray': False,
	'mathtext.bf': 'serif:bold',
	'mathtext.cal': 'cursive',
	'mathtext.fallback_to_cm': True,
	'mathtext.fontset': 'cm',
	'mathtext.it': 'serif:italic',
	'mathtext.rm': 'serif',
	'mathtext.sf': 'sans\\-serif',
	'mathtext.tt': 'monospace',
	'numerix': 'numpy',
	'patch.antialiased': True,
	'patch.edgecolor': 'k',
	'patch.facecolor': 'b',
	'patch.linewidth': 1.0,
	'path.simplify': False,
	'pdf.compression': 6,
	'pdf.fonttype': 3,
	'pdf.inheritcolor': False,
	'pdf.use14corefonts': False,
	'plugins.directory': '.matplotlib_plugins',
	'polaraxes.grid': True,
	'ps.distiller.res': 6000,
	'ps.fonttype': 3,
	'ps.papersize': 'letter',
	'ps.useafm': False,
	'ps.usedistiller': False,
	'savefig.dpi': 100,
	'savefig.edgecolor': 'w',
	'savefig.facecolor': 'w',
	'savefig.orientation': 'portrait',
	'svg.embed_char_paths': True,
	'svg.image_inline': True,
	'svg.image_noscale': False,
	'text.color': 'k',
	'text.dvipnghack': None,
	'text.fontangle': 'normal',
	'text.fontsize': 'x-large',
	'text.fontstyle': 'normal',
	'text.fontvariant': 'normal',
	'text.fontweight': 'normal',
	'text.latex.preamble': [''],
	'text.latex.unicode': False,
	'text.usetex': False,
	'timezone': 'UTC',
	'tk.pythoninspect': False,
	'tk.window_focus': False,
	'toolbar': 'toolbar2',
	'units': False,
	'verbose.fileo': 'sys.stdout',
	'verbose.level': 'silent',
	'xtick.color': 'k',
	'xtick.direction': 'out', # 'in'
	'xtick.labelsize': 'x-large', # 'medium'
	'xtick.major.pad': 4,
	'xtick.major.size': 4,
	'xtick.minor.pad': 4,
	'xtick.minor.size': 2,
	'ytick.color': 'k',
	'ytick.direction': 'out', # 'in'
	'ytick.labelsize': 'x-large', # 'medium'
	'ytick.major.pad': 4,
	'ytick.major.size': 4,
	'ytick.minor.pad': 4,
	'ytick.minor.size': 2
}



def add_subplot(fig = None, *args):
	if fig != None:
		figure(fig)
	gcf().add_subplot(*args)
	draw()


def get_axis(fig = None, N = None):
	if fig != None:
		figure(fig)
	"""
	function to handle axis in complexe subplot.
	if N == None, return current axis, eq. gca()
	if N != None, return selected (or last axis)
	"""
	
	if N == None:
		print("Axes length: " + str(len(gcf().axes)))
		return gcf().axes[-1]
	else:
		if N > len(gcf().axes) - 1:
			print("Index out of bounds")
			return gcf().axes[-1]
		else:
			return gcf().axes[N]

class graph:
	"""
	Overlaod of figure() + plot() function
	a 2 in 1
	graph(...) to open and plot
	plot(...) to append a line in the graph
	"""
	
	f = None
	p = None
	def __init__(self, *args):
		self.f = figure()
		# add a fancy color on the graph's border
		# self.f.set_facecolor("#dffa87")
		self.p = plot(*args)


class add_plot:
	def __init__(self, fig = None, *args):
		if fig != None:
			self.f = figure(fig)
		self.p = plot(*args)

	


class imgraph():
	"""
	Overlaod of figure() + imshow() function
	a 2 in 1
	imgraph(...) to open and plot
	imshow(...) to append a line in the graph
	"""

	f = None
	i = None
	def __init__(self, X, cmap=None, norm=None, aspect=None, interpolation=None, alpha=1.0, vmin=None, vmax=None, origin="lower", extent=None, **kwargs):
		self.f = figure()
		# another fancy color
		# self.f.set_facecolor("#fcd628")
		self.i = imshow(X, cmap=cmap, norm=norm, aspect=aspect, interpolation=interpolation, alpha=alpha, vmin=vmin, vmax=vmax, origin=origin, extent=extent, **kwargs)



class colorgraph():
	"""
	Overlaod of figure() + pcolor() function
	a 2 in 1
	colorgraph(...) to open and plot
	pcolor(...) to append a line in the graph
	"""

	f = None
	i = None
	def __init__(self, *args):
		self.f = figure()
		#self.f.set_facecolor("#fcd628")
		self.i = pcolormesh(*args)



class JSONgraph():
	"""
	Whole figure in JSON.
	"""

	figure = None
	figure_dict = {}

	def __init__(self, fig = None):
		if fig != None:
			figure(fig)
		self.figure = gcf()
		
		self.figure_dict = {}
		
		self.figure_dict["alpha"] = self.figure.get_alpha()
		self.figure_dict["animated"] = self.figure.get_animated()
		self.figure_dict["clip_on"] = self.figure.get_clip_on()
		self.figure_dict["dpi"] = self.figure.get_dpi()
		self.figure_dict["edgecolor"] = self.figure.get_edgecolor()
		self.figure_dict["facecolor"] = self.figure.get_facecolor()
		self.figure_dict["figheight"] = self.figure.get_figheight()
		self.figure_dict["figwidth"] = self.figure.get_figwidth()
		self.figure_dict["frameon"] = self.figure.get_frameon()
		self.figure_dict["label"] = self.figure.get_label()
		self.figure_dict["visible"] = self.figure.get_visible()

		self.figure_dict["axes"] = []

		if len(self.figure.axes) > 0:
			for _axis in self.figure.axes:

				_axis_dict = {}

				_axis_dict["adjustable"] = _axis.get_adjustable()
				_axis_dict["alpha"] = _axis.get_alpha()
				_axis_dict["anchor"] = _axis.get_anchor()
				_axis_dict["animated"] = _axis.get_animated()
				_axis_dict["aspect"] = _axis.get_aspect()
				_axis_dict["autoscale_on"] = _axis.get_autoscale_on()
				_axis_dict["axis_bgcolor"] = _axis.get_axis_bgcolor()
				_axis_dict["axisbelow"] = _axis.get_axisbelow()
				_axis_dict["clip_on"] = _axis.get_clip_on()
				_axis_dict["frame_on"] = _axis.get_frame_on()
				_axis_dict["label"] = _axis.get_label()
				_axis_dict["position"] = _axis.get_position().bounds
				_axis_dict["title"] = _axis.get_title()
				_axis_dict["visible"] = _axis.get_visible()
				_axis_dict["xlabel"] = _axis.get_xlabel()
				_axis_dict["xscale"] = _axis.get_xscale()
				_axis_dict["ylabel"] = _axis.get_ylabel()
				_axis_dict["yscale"] = _axis.get_yscale()
				_axis_dict["zorder"] = _axis.get_zorder()

				_axis_dict["lines"] = []
				_axis_dict["texts"] = []

				if len(_axis.lines) > 0:
					for _line in _axis.lines:
						
						_line_dict = {}

						_line_dict["alpha"] = _line.get_alpha()
						_line_dict["animated"] = _line.get_animated()
						_line_dict["antialiased"] = _line.get_antialiased()
						_line_dict["clip_on"] = _line.get_clip_on()
						_line_dict["color"] = _line.get_color()
						_line_dict["dash_capstyle"] = _line.get_dash_capstyle()
						_line_dict["dash_joinstyle"] = _line.get_dash_joinstyle()
						_line_dict["drawstyle"] = _line.get_drawstyle()
						_line_dict["label"] = _line.get_label()
						_line_dict["linestyle"] = _line.get_linestyle()
						_line_dict["linewidth"] = _line.get_linewidth()
						_line_dict["marker"] = _line.get_marker()
						_line_dict["markeredgecolor"] = _line.get_markeredgecolor()
						_line_dict["markeredgewidth"] = _line.get_markeredgewidth()
						_line_dict["markerfacecolor"] = _line.get_markerfacecolor()
						_line_dict["markersize"] = _line.get_markersize()
						_line_dict["solid_capstyle"] = _line.get_solid_capstyle()
						_line_dict["solid_joinstyle"] = _line.get_solid_joinstyle()
						_line_dict["visible"] = _line.get_visible()
						_line_dict["xdata"] = list(_line.get_xdata())
						_line_dict["ydata"] = list(_line.get_ydata())
						_line_dict["zorder"] = _line.get_zorder()

						_axis_dict["lines"].append(_line_dict)
				
				if len(_axis.texts) > 0:
					for _text in _axis.texts:

						_text_dict = {}

						_text_dict["alpha"] = _text.get_alpha()
						_text_dict["animated"] = _text.get_animated()
						_text_dict["bbox_patch"] = _text.get_bbox_patch()
						_text_dict["clip_on"] = _text.get_clip_on()
						_text_dict["color"] = _text.get_color()
						_text_dict["family"] = _text.get_family()
						_text_dict["figure"] = _text.get_figure()
						_text_dict["horizontalalignment"] = _text.get_horizontalalignment()
						_text_dict["label"] = _text.get_label()
						_text_dict["name"] = _text.get_name()
						_text_dict["position"] = _text.get_position()
						_text_dict["rotation"] = _text.get_rotation()
						_text_dict["size"] = _text.get_size()
						_text_dict["stretch"] = _text.get_stretch()
						_text_dict["style"] = _text.get_style()
						_text_dict["text"] = _text.get_text()
						_text_dict["variant"] = _text.get_variant()
						_text_dict["verticalalignment"] = _text.get_verticalalignment()
						_text_dict["visible"] = _text.get_visible()
						_text_dict["weight"] = _text.get_weight()
						_text_dict["zorder"] = _text.get_zorder()

						_axis_dict["texts"].append(_text_dict)

				self.figure_dict["axes"].append(_axis_dict)


	
	def save(self, fname = None):
		if fname == None:
			print("Filename needed")
			return
		else:
			_f = file(fname, "w")
			sj.dump(self.figure_dict, _f)
			_f.close()




class openNanoQtgraph():
    g = None
    ax = None
    def __init__(self, dictdata, reversed = False):
        self.g = graph()
        draw_legend = False;
        self.ax = gca()
        
        self.ax.set_title(dictdata["title"])
        self.ax.set_xlabel(dictdata["x_label"])
        self.ax.set_ylabel(dictdata["y_label"])
    
        if dictdata["logscale_x"]:
            self.ax.set_xscale("log")
        else:
            self.ax.set_xscale("linear")
        
        if dictdata["logscale_y"]:
            self.ax.set_yscale("log")
        else:
            self.ax.set_yscale("linear")
    
        if dictdata.has_key("curves"):
            for n in range(len(dictdata["curves"])):
                if dictdata["curves"][n] != None:
                    curvex, curvey = hsplit(array(dictdata["curves"][n]["data"]), 2)
                    plot(curvex, curvey)
                    curveline = self.ax.lines[-1]
                    if dictdata["curves"][n]["options"].has_key("label"):
                        curveline.set_label(dictdata["curves"][n]["options"]["label"])
                        draw_legend = True; 
                    curveline.set_linestyle(curve_line[dictdata["curves"][n]["options"]["pen_style"]])
                    curveline.set_linewidth(dictdata["curves"][n]["options"]["pen_width"])
                    if str(dictdata["curves"][n]["options"]["pen_color"]).startswith("0x"):
                        tmp_col = hex_x_to_normalized_tuple(str(dictdata["curves"][n]["options"]["pen_color"]))
                        curveline.set_color(tmp_col)# color_list[n%7])
                        curveline.set_markeredgecolor(tmp_col)# color_list[n%7])
                        curveline.set_markerfacecolor(tmp_col)# color_list[n%7])
                    else:
                        curveline.set_color(dictdata["curves"][n]["options"]["pen_color"])
                        curveline.set_markeredgecolor(dictdata["curves"][n]["options"]["pen_color"])
                        curveline.set_markerfacecolor(dictdata["curves"][n]["options"]["pen_color"])

                    curveline.set_marker(curve_marker[dictdata["curves"][n]["options"]["symbol"]])
                    curveline.set_markersize(dictdata["curves"][n]["options"]["symbol_size"])
        
        if dictdata.has_key("array"):
            a = array(dictdata["array"]["data"]) 
            x, y = meshgrid(linspace(dictdata["x_min"], dictdata["x_max"], len(a)), linspace(dictdata["y_min"], dictdata["y_max"], len(a[0])))
            if reversed:
                pcolormesh(x, y, a.max() - a)
            else:
                pcolormesh(x, y, a)
        
        if draw_legend:
            legend() 

        draw()


class openNanoQtgraph_polar():
    g = None
    ax = None
    def __init__(self, dictdata, reversed = False, deg_rad = 1):
        figure()    
        polar()
        if dictdata.has_key("array"):
            a = array(dictdata["array"]["data"]) 
            x, y = meshgrid(linspace(dictdata["x_min"], dictdata["x_max"], len(a)), linspace(dictdata["y_min"], dictdata["y_max"], len(a[0])))
            if reversed:
                pcolormesh(deg_rad*y, x, a.max() - a)
            else:
                pcolormesh(deg_rad*y, x, a)

        draw()


class openTSVgraph():
	def __init__(self, fname):
		_data = csv2rec(fname, delimiter=" ")
		graph()
		for i in _data.dtype.names:
			plot(_data[i])	
			gca().lines[-1].set_label(i)
		legend()
		title(fname)
		xlabel("i")
		ylabel("var")


################################################################################

class SetupAxes(eta.HasTraits):
	
    axes_list_max = eta.Str
    axes_list_index = eta.Int(-1)
    get_current_axes = eta.Button("Get Current Axes")
    axis = None
    xbound = eta.Str
    xlabel = eta.Str
    xfontsize = eta.Enum('large', 'xx-small', 'x-small', 'small', 'medium', 'x-large', 'xx-large')
    xfontstyle = eta.Enum('normal', 'italic', 'oblique')
    xfontweight = eta.Enum('normal', 'ultralight', 'light', 'regular', 'book', 'medium', 'roman', 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black')
    xscale = eta.Enum("linear", "log")
    ybound = eta.Str
    ylabel = eta.Str
    yfontsize = eta.Enum('large', 'xx-small', 'x-small', 'small', 'medium', 'x-large', 'xx-large')
    yfontstyle = eta.Enum('normal', 'italic', 'oblique')
    yfontweight = eta.Enum('normal', 'ultralight', 'light', 'regular', 'book', 'medium', 'roman', 'semibold', 'demibold', 'demi', 'bold', 'heavy', 'extra bold', 'black')
    yscale = eta.Enum("linear", "log")
    adjustable = eta.Enum("box", "datalim")
    alpha = eta.Range(0.0, 1.0) # RANGE
    anchor = eta.Enum("C", "SW", "S", "SE", "E", "NE", "N", "NW", "W")
    animated = eta.Bool # BOOL
    aspect = eta.Enum("auto", "normal", "equal")
    autoscale_on = eta.Bool # BOOL
    axes = eta.Str
    axisbelow = eta.Bool # BOOL
    axis_bgcolor = eta.Color
    clip_box = eta.Str
    clip_on = eta.Bool # BOOL
    clip_path = eta.Str
    contains = eta.Str
    cursor_props = eta.Str
    figure = eta.Str
    frame_on = eta.Bool # BOOL
    label = eta.Str
    navigate = eta.Bool # BOOL
    navigate_mode = eta.Str
    picker = eta.Bool # BOOL
    position_left = eta.Range(0.0, 1.0)
    position_bottom = eta.Range(0.0, 1.0)
    position_width = eta.Range(0.0, 1.0)
    position_height = eta.Range(0.0, 1.0)
    rasterization_zorder = eta.Str
    rasterized = eta.Str
    snap = eta.Str
    title = eta.Str
    transform = eta.Str
    url = eta.Str
    visible = eta.Bool # BOOL
    zorder = eta.Float

    def _update(self):
        self.xbound = strlist(self.axis.get_xbound())
        self.xlabel = self.axis.get_xlabel()
        self.xscale = self.axis.get_xscale()
        self.ybound = strlist(self.axis.get_ybound())
        self.ylabel = self.axis.get_ylabel()
        self.yscale = self.axis.get_yscale()
        self.adjustable = self.axis.get_adjustable()
        self.alpha = self.axis.get_alpha()
        self.anchor = self.axis.get_anchor()
        self.animated = self.axis.get_animated()
        self.aspect = self.axis.get_aspect()
        self.autoscale_on = self.axis.get_autoscale_on()
        self.axes = str(self.axis.get_axes())
        self.axisbelow = self.axis.get_axisbelow()
	if type(self.axis.get_axis_bgcolor()) == tuple:
		self.axis_bgcolor = (
			int(255*self.axis.get_axis_bgcolor()[0]),
			int(255*self.axis.get_axis_bgcolor()[1]),
			int(255*self.axis.get_axis_bgcolor()[2])
			)
	elif self.axis.get_axis_bgcolor() in color_corr.keys():
		self.axis_bgcolor = color_corr[self.axis.get_axis_bgcolor()]
	else:
		self.axis_bgcolor = self.axis.get_axis_bgcolor()
        self.clip_box = str(self.axis.get_clip_box())
        self.clip_on = self.axis.get_clip_on()
        self.clip_path = str(self.axis.get_clip_path())
        self.contains = str(self.axis.get_contains())
        self.cursor_props = strlist(self.axis.get_cursor_props())
        self.figure = str(self.axis.get_figure())
        self.frame_on = self.axis.get_frame_on()
        self.label = self.axis.get_label()
        self.navigate = self.axis.get_navigate()
        self.navigate_mode = str(self.axis.get_navigate_mode())
        if self.axis.get_picker():
            self.picker = False
        else:
            self.picker = True
	tmp_position = self.axis.get_position().bounds
        self.position_left = tmp_position[0]
        self.position_bottom = tmp_position[1]
        self.position_width = tmp_position[2]
        self.position_height = tmp_position[3]
        self.snap = str(self.axis.get_snap())
        self.title = self.axis.get_title()
        self.transform = str(self.axis.get_transform())
        self.url = str(self.axis.get_url())
        self.visible = self.axis.get_visible()
	self.zorder = self.axis.get_zorder()

    def __init__(self):
        self.axes_list_max = strlist(range(len(gcf().axes)))
	self.axis = gcf().axes[self.axes_list_index]
        self._update()

    def _get_current_axes_fired(self):
        self.axes_list_max = strlist(range(len(gcf().axes)))
	self.axis = gcf().axes[self.axes_list_index]
        self._update()
	
    def _axes_list_index_changed(self):
	self.axis = gcf().axes[self.axes_list_index]
        self._update()

        
    def _alpha_changed(self):
        self.axis.set_alpha(self.alpha)
        draw()

    def _aspect_changed(self):
        self.axis.set_aspect(self.aspect)
        draw()
            
    def _autoscale_on_changed(self):
        self.axis.set_autoscale_on(self.autoscale_on)
        draw()
            
    def _axisbelow_changed(self):
        self.axis.set_axisbelow(self.axisbelow)
        draw()
            
    def _axis_bgcolor_changed(self):
        self.axis.set_axis_bgcolor(
		(
			self.axis_bgcolor[0]/255.,
			self.axis_bgcolor[1]/255.,
			self.axis_bgcolor[2]/255.
			)
		)
        draw()
            
    def _clip_on_changed(self):
        self.axis.set_clip_on(self.clip_on)
        draw()
            
    def _label_changed(self):
        self.axis.set_label(self.label)
        draw()
            
    def _position_left_changed(self):
        self.axis.set_position([self.position_left, self.position_bottom, self.position_width, self.position_height])
        draw()

    def _position_bottom_changed(self):
        self.axis.set_position([self.position_left, self.position_bottom, self.position_width, self.position_height])
        draw()

    def _position_width_changed(self):
        self.axis.set_position([self.position_left, self.position_bottom, self.position_width, self.position_height])
        draw()

    def _position_height_changed(self):
        self.axis.set_position([self.position_left, self.position_bottom, self.position_width, self.position_height])
        draw()
            
    def _title_changed(self):
        self.axis.set_title(self.title)
        draw()
            
    def _visible_changed(self):
        self.axis.set_visible(self.visible)
        draw()
            
    def _xlabel_changed(self):
        self.axis.set_xlabel(self.xlabel, size=self.xfontsize, style=self.xfontstyle, weight=self.xfontweight)
        draw()

    def _ylabel_changed(self):
        self.axis.set_ylabel(self.ylabel, size=self.yfontsize, style=self.yfontstyle, weight=self.yfontweight)
        draw()
            
    def _xscale_changed(self):
        self.axis.set_xscale(self.xscale)
        draw()

    def _yscale_changed(self):
        self.axis.set_yscale(self.yscale)
        draw()

    def _adjustable_changed(self):
        self.axis.set_adjustable(self.adjustable)
        draw()

    def _anchor_changed(self):
        self.axis.set_anchor(self.anchor)
        draw()

    def _animated_changed(self):
        self.axis.set_animated(self.animated)
        draw()

    def _frame_on_changed(self):
        self.axis.set_frame_on(self.frame_on)
        draw()

    def _navigate_changed(self):
        self.axis.set_navigate(self.navigate)
        draw()

    def _picker_changed(self):
        self.axis.set_picker(self.picker)
        draw()

    def _zorder_changed(self):
        self.axis.set_zorder(self.zorder)
        draw()

        
    traits_view = etua.View(
	etua.Item("get_current_axes", show_label=False),
	etua.Item("axes_list_max", style="readonly"),
        etua.Item("axes_list_index"),
        etua.Item("axes", width=640),
        etua.Group(
            etua.Group(
		"title",
                "alpha",
                "aspect",
                "frame_on",
                "visible",
                "autoscale_on",
		"axisbelow",
                "clip_on",
                etua.Item("position_left"),
                etua.Item("position_bottom"),
                etua.Item("position_width"),
                etua.Item("position_height"),
		"zorder",
                etua.Item("axis_bgcolor", style="custom"),
                label="general", dock="tab"),
            etua.Group(
                "xlabel",
                "xfontsize",
                "xfontstyle",
                "xfontweight",
                "xbound",
                "xscale",
                "_",
                "ylabel",
                "yfontsize",
                "yfontstyle",
                "yfontweight",
                "ybound",
                "yscale",
                label="X Y", dock="tab"),
            etua.Group(
                "adjustable",
                "anchor",
                "animated",
                "axes",
                "axes_locator",
                "clip_box",
                "clip_path",
                "contains",
                "cursor_props",
                "figure",
                "gid",
                "label",
                "navigate",
                "navigate_mode",
                "picker",
		"snap",
                "transform",
                "url",
                label="other", dock="tab"),
            layout="tabbed"
            ),
        resizable = True,
	scrollable = True
        )


def setaxes():
	SetupAxes().configure_traits()



#######################################################################################

class SetupLines(eta.HasTraits):
	line_list_max = eta.Str
	line_list_index = eta.Int(-1)
	get_current_lines = eta.Button("Get Current Lines")
	line = None
	show_all = eta.Bool(True)
	alpha = eta.Range(0.0, 1.0)
	animated = eta.Bool
	antialiased = eta.Bool
	axes = eta.Str
	children = eta.Str
	clip_box = eta.Str
	clip_on = eta.Bool
	clip_path = eta.Str
	color = eta.Color
	contains = eta.Str
	dash_capstyle = eta.Enum("butt", "round", "projecting")
	dash_joinstyle = eta.Enum("miter", "round", "bevel")
	data = eta.Str
	drawstyle = eta.Enum("default", "steps", "steps-pre", "steps-mid", "steps-post")
	figure = eta.Str
	fillstyle = eta.Enum("full", "left", "right", "bottom", "top")
	gid = eta.Str
	label = eta.Str
	linestyle = eta.Enum("None", "-", "--", "-.", ":")
	linewidth = eta.Range(0., 100.)
	marker = eta.Enum("None", ".", ",", "o", "v", "^", "<", ">", "1", "2", "3", "4", "s", "p", "*", "h", "H", "+", "x", "D", "d", "_")
	markeredgecolor = eta.Color
	markeredgewidth = eta.Range(0., 100.)
	markerfacecolor = eta.Color
	markersize = eta.Range(0., 100.0)
	markevery = eta.Int
	path = eta.Str
	picker = eta.Float
	pickradius = eta.Range(0.1, 100.0)
	rasterized = eta.Str
	snap = eta.Str
	solid_capstyle = eta.Enum("butt", "round", "projecting")
	solid_joinstyle = eta.Enum("miter", "round", "bevel")
	transform = eta.Str
	transformed_clip_path_and_affine = eta.Str
	url = eta.Str
	visible = eta.Bool
	xdata = eta.Str
	xydata = eta.Str
	ydata = eta.Str
	zorder = eta.Int

	def _update(self):
		if len(gca().lines) == 0:
			plot(rand(10), rand(10))
		self.line_list_max = strlist(range(len(gca().lines)))
		if abs(self.line_list_index) > len(gca().lines) - 1:
			self.line_list_index = len(gca().lines) - 1
		self.line = gca().lines[self.line_list_index]
		self.alpha = self.line.get_alpha()
		self.animated = self.line.get_animated()
		self.antialiased = self.line.get_antialiased()
		self.axes = str(self.line.get_axes())
		self.children = strlist(self.line.get_children())
		self.clip_box = str(self.line.get_clip_box())
		self.clip_on = self.line.get_clip_on()
		self.clip_path = str(self.line.get_clip_path())
		if type(self.line.get_color()) == tuple:
			self.color = (
				int(255*self.line.get_color()[0]),
				int(255*self.line.get_color()[1]),
				int(255*self.line.get_color()[2])
				)
		elif self.line.get_color() in color_corr.keys():
			self.color = color_corr[self.line.get_color()]
		else:
			self.color = self.line.get_color()
		self.contains = str(self.line.get_contains())
		self.dash_capstyle = self.line.get_dash_capstyle()
		self.dash_joinstyle = self.line.get_dash_joinstyle()
		self.data = strlist(self.line.get_data())
		self.drawstyle = self.line.get_drawstyle()
		self.figure = str(self.line.get_figure())
		self.label = self.line.get_label()
		self.linestyle = self.line.get_linestyle()
		self.linewidth = self.line.get_linewidth()
		self.marker = self.line.get_marker()
		if type(self.line.get_markeredgecolor()) == tuple:
			self.markeredgecolor = (
				int(255*self.line.get_markeredgecolor()[0]),
				int(255*self.line.get_markeredgecolor()[1]),
				int(255*self.line.get_markeredgecolor()[2])
				)
		elif self.line.get_markeredgecolor() in color_corr.keys():
			self.markeredgecolor = color_corr[self.line.get_markeredgecolor()]
		else:
			self.markeredgecolor = self.line.get_markeredgecolor()
			
		if type(self.line.get_markerfacecolor()) == tuple:
			self.markerfacecolor = (
				int(255*self.line.get_markerfacecolor()[0]),
				int(255*self.line.get_markerfacecolor()[1]),
				int(255*self.line.get_markerfacecolor()[2])
				)
		elif self.line.get_markerfacecolor() in color_corr.keys():
			self.markerfacecolor = color_corr[self.line.get_markerfacecolor()]
		else:
			self.markerfacecolor = self.line.get_markerfacecolor()
		self.markeredgewidth = self.line.get_markeredgewidth()
		self.markersize = self.line.get_markersize()
		self.markevery = 1
		self.path = str(self.line.get_path())
		if self.line.get_picker() == None:
			self.picker = False
		else:
			self.picker = self.line.get_picker()
		self.snap = str(self.line.get_snap())
		self.solid_capstyle = self.line.get_solid_capstyle()
		self.solid_joinstyle = self.line.get_solid_joinstyle()
		self.transform = str(self.line.get_transform())
		self.transformed_clip_path_and_affine = str(self.line.get_transformed_clip_path_and_affine())
		self.url = str(self.line.get_url())
		self.visible = self.line.get_visible()
		self.xdata = strlist(self.line.get_xdata())
		self.xydata = strlist(self.line.get_xydata())
		self.ydata = strlist(self.line.get_ydata())
		self.zorder = self.line.get_zorder()

	def __init__(self):
		self._update()

	def _get_current_lines_fired(self):
		self.show_all = True
		if len(gca().lines) == 0:
			plot(rand(10), rand(10))
		self.line_list_max = strlist(range(len(gca().lines)))
		if abs(self.line_list_index) > len(gca().lines) - 1:
			self.line_list_index = len(gca().lines) - 1
		self.line = gca().lines[self.line_list_index]
		self._update()
		draw()

	def _line_list_index_changed(self):
		self.line_list_max = strlist(range(len(gca().lines)))
		if abs(self.line_list_index) > len(gca().lines) - 1:
			self.line_list_index = len(gca().lines) - 1
		self.line = gca().lines[self.line_list_index]
		if self.show_all == False:
			for i in range(len(gca().lines)):
				gca().lines[i].set_visible(False)
			gca().lines[self.line_list_index].set_visible(True)
		else:
			for i in range(len(gca().lines)):
				gca().lines[i].set_visible(True)
		self._update()
		draw()

	def _show_all_changed(self):
		if self.show_all == False:
			for i in range(len(gca().lines)):
				gca().lines[i].set_visible(False)
			gca().lines[self.line_list_index].set_visible(True)
		else:
			for i in range(len(gca().lines)):
				gca().lines[i].set_visible(True)
		draw()

	def _alpha_changed(self):
		self.line.set_alpha(self.alpha)
		draw()

	def _animated_changed(self):
		self.line.set_animated(self.animated)
		draw()

	def _antialiased_changed(self):
		self.line.set_antialiased(self.antialiased)
		draw()

	def _clip_on_changed(self):
		self.line.set_clip_on(self.clip_on)
		draw()

	def _color_changed(self):
		self.line.set_color((self.color[0]/255., self.color[1]/255., self.color[2]/255.))
		draw()

	def _dash_capstyle_changed(self):
		self.line.set_dash_capstyle(self.dash_capstyle)
		draw()

	def _dash_joinstyle_changed(self):
		self.line.set_dash_joinstyle(self.dash_joinstyle)
		draw()

	def _drawstyle_changed(self):
		self.line.set_drawstyle(self.drawstyle)
		draw()

	def _label_changed(self):
		self.line.set_label(self.label)
		draw()

	def _linestyle_changed(self):
		self.line.set_linestyle(self.linestyle)
		draw()

	def _linewidth_changed(self):
		self.line.set_linewidth(self.linewidth)
		draw()

	def _marker_changed(self):
		self.line.set_marker(self.marker)
		draw()

	def _markeredgecolor_changed(self):
		self.line.set_markeredgecolor(
			(
				self.markeredgecolor[0]/255.,
				self.markeredgecolor[1]/255.,
				self.markeredgecolor[2]/255.
				)
			)
		draw()

	def _markeredgewidth_changed(self):
		self.line.set_markeredgewidth(self.markeredgewidth)
		draw()

	def _markerfacecolor_changed(self):
		self.line.set_markerfacecolor(
			(
				self.markerfacecolor[0]/255.,
				self.markerfacecolor[1]/255.,
				self.markerfacecolor[2]/255.
				)
			)
		draw()

	def _markersize_changed(self):
		self.line.set_markersize(self.markersize)
		draw()

	def _picker_changed(self):
		self.line.set_picker(self.picker)
		draw()

	def _snap_changed(self):
		self.line.set_snap(self.snap)
		draw()

	def _solid_capstyle_changed(self):
		self.line.set_solid_capstyle(self.solid_capstyle)
		draw()

	def _solid_joinstyle_changed(self):
		self.line.set_solid_joinstyle(self.solid_joinstyle)
		draw()

	def _url_changed(self):
		self.line.set_url(self.url)
		draw()

	def _visible_changed(self):
		self.line.set_visible(self.visible)
		draw()

	def _xdata_changed(self):
		self.line.set_xdata(eval(self.xdata))
		draw()

	def _ydata_changed(self):
		self.line.set_ydata(eval(self.ydata))
		draw()

	def _zorder_changed(self):
		self.line.set_zorder(self.zorder)
		draw()


	traits_view = etua.View(
		etua.Item("get_current_lines", show_label = False),
		etua.Item("line_list_max", style="readonly"),
		"line_list_index",
		"show_all",
		etua.Item("label", width=640),
		etua.Group(
			etua.Group(
				etua.Item("color"), #style="custom"),
				"dash_capstyle",
				"dash_joinstyle",
				"solid_capstyle",
				"solid_joinstyle",
				"drawstyle",
				"linestyle",
				"linewidth",
				"_",
				etua.Item("markeredgecolor"), #style="custom"),
				etua.Item("markerfacecolor"), #style="custom"),
				"marker",
				"markeredgewidth",
				"markersize",
				label="line & marker", dock="tab"),
			etua.Group(
				"alpha",
				"antialiased",
				"clip_on",
				"label",
				"picker",
				"visible",
				"zorder",
				label="general", dock="tab"),
			etua.Group(
				"data",
				"xdata",
				"ydata",
				"xydata",
				label="data", dock="tab"),
			etua.Group(
				"animated",
				"axes",
				"children",
				"clip_box",
				"clip_path",
				"contains",
				"figure",
				"path",
				"snap",
				"transform",
				"transformed_clip_path_and_affine",
				"url",
				"window_extent",
				label="other", dock="tab"),
			layout="tabbed"),
		resizable = True,
		scrollable = True
		)
    
def setlines():
        SetupLines().configure_traits()






###############################################################################

class SetupText(eta.HasTraits):
	t_list_max = eta.Any
	t_list_index = eta.Int(-1)
	t = None
	get_current_texts = eta.Button("Get Current Texts")
	alpha = eta.Range(0.0, 1.0)
	animated = eta.Bool
	bbox_facecolor = eta.Color
	bbox_alpha = eta.Range(0.0, 1.0)
	clip_box = eta.Any
	clip_on = eta.Bool
	clip_path = eta.Any
	color = eta.Color
	family = eta.Any(desc="FONTNAME, 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'")
	horizontalalignment = eta.Enum('center', 'right', 'left')
	label = eta.Any
	name = eta.Any(desc="FONTNAME, 'serif', 'sans-serif', 'cursive', 'fantasy', 'monospace'")
	picker = eta.Any
	position_x = eta.Float
	position_y = eta.Float
	#prop_tup = eta.Any
	rotation = eta.Str("0", desc="angle in degrees | 'vertical' | 'horizontal'")
	size = eta.Any(desc="size in points | 'xx-small' | 'x-small' | 'small' | 'medium' | 'large' | 'x-large' | 'xx-large'")
	snap = eta.Any
	stretch = eta.Any(desc="a numeric value in range 0-1000 | 'ultra-condensed' | 'extra-condensed' | 'condensed' | 'semi-condensed' | 'normal' | 'semi-expanded' | 'expanded' | 'extra-expanded' | 'ultra-expanded'")
	style = eta.Enum('normal', 'italic', 'oblique')
	text = eta.Any
	url = eta.Any
	variant = eta.Enum('normal', 'small-caps')
	verticalalignment = eta.Enum('center', 'top', 'bottom', 'baseline')
	visible = eta.Bool
	weight = eta.Any(desc="a numeric value in range 0-1000 | 'ultralight' | 'light' | 'normal' | 'regular' | 'book' | 'medium' | 'roman' | 'semibold' | 'demibold' | 'demi' | 'bold' | 'heavy' | 'extra bold' | 'black'")
	zorder = eta.Float

	def __init__(self):
		self.t_list_max = strlist(range(len(gca().texts)))
		self.t_list_index = -1
		self.t = gca().texts[self.t_list_index]
		for i in gca().texts:
			if i._bbox == None:
				i._bbox = {"facecolor": "white", "alpha": 1.0}
		self._update()

	def _get_current_texts_fired(self):
		self.t_list_max = strlist(range(len(gca().texts)))
		self.t_list_index = -1
		self.t = gca().texts[self.t_list_index]
		for i in gca().texts:
			if i._bbox == None:
				i._bbox = {"facecolor": "white", "alpha": 1.0}
		self._update()

	def _t_list_index_changed(self):
		self.t = gca().texts[self.t_list_index]
		self._update()

	def _update(self):
		self.alpha = self.t.get_alpha()
		self.animated = self.t.get_animated()
		if type(self.t._bbox["facecolor"]) == tuple:
			self.bbox_facecolor = (
				int(255*self.t._bbox["facecolor"][0]),
				int(255*self.t._bbox["facecolor"][1]),
				int(255*self.t._bbox["facecolor"][2])
				)
		elif self.t._bbox["facecolor"] in color_corr.keys():
			self.bbox_facecolor = color_corr[self.t._bbox["facecolor"]]
		else:
			self.bbox_facecolor = self.t._bbox["facecolor"]
		self.bbox_alpha = self.t._bbox["alpha"]
		self.clip_box = self.t.get_clip_box()
		self.clip_on = self.t.get_clip_on()
		self.clip_path = self.t.get_clip_path()
		if type(self.t.get_color()) == tuple:
			self.color = (
				int(255*self.t.get_color()[0]),
				int(255*self.t.get_color()[1]),
				int(255*self.t.get_color()[2])
				)
		elif self.t.get_color() in color_corr.keys():
			self.color = color_corr[self.t.get_color()]
		else:
			self.color = self.t.get_color()
		self.family = self.t.get_family()
		self.horizontalalignment = self.t.get_horizontalalignment()
		self.label = self.t.get_label()
		self.name = self.t.get_name()
		self.picker = self.t.get_picker()
		tmp_pos = self.t.get_position()
		self.position_x = tmp_pos[0]
		self.position_y = tmp_pos[1]
		#self.prop_tup = self.t.get_prop_tup()
		self.rotation = str(self.t.get_rotation())
		self.size = self.t.get_size()
		self.snap = self.t.get_snap()
		self.stretch = self.t.get_stretch()
		self.style = self.t.get_style()
		self.text = self.t.get_text()
		self.url = self.t.get_url()
		self.variant = self.t.get_variant()
		self.verticalalignment = self.t.get_verticalalignment()
		self.visible = self.t.get_visible()
		self.weight = self.t.get_weight()
		self.zorder = float(self.t.get_zorder())

	def _alpha_changed(self):
		self.t.set_alpha(self.alpha)
		draw()

	def _animated_changed(self):
		self.t.set_animated(self.animated)
		draw()

	def _bbox_facecolor_changed(self):
		self.t._bbox.update(
			{
				"facecolor": (
					self.bbox_facecolor[0]/255.,
					self.bbox_facecolor[1]/255.,
					self.bbox_facecolor[2]/255.
					)
				}
			)
		draw()

	def _bbox_alpha_changed(self):
		self.t._bbox.update({"alpha": self.bbox_alpha})
		draw()

	def _clip_box_changed(self):
		self.t.set_clip_box(self.clip_box)
		draw()

	def _clip_on_changed(self):
		self.t.set_clip_on(self.clip_on)
		draw()

	def _clip_path_changed(self):
		self.t.set_clip_path(self.clip_path)
		draw()

	def _color_changed(self):
		self.t.set_color(
			(
				self.color[0]/255.,
				self.color[1]/255.,
				self.color[2]/255.
				)
			)
		draw()

	def _family_changed(self):
		self.t.set_family(self.family)
		draw()

	def _horizontalalignment_changed(self):
		self.t.set_horizontalalignment(self.horizontalalignment)
		draw()

	def _label_changed(self):
		self.t.set_label(self.label)
		draw()

	def _name_changed(self):
		self.t.set_name(self.name)
		draw()

	def _picker_changed(self):
		self.t.set_picker(self.picker)
		draw()

	def _position_x_changed(self):
		self.t.set_x(self.position_x)
		draw()

	def _position_y_changed(self):
		self.t.set_y(self.position_y)
		draw()

	def _rotation_changed(self):
		self.t.set_rotation(eval(self.rotation))
		draw()

	def _size_changed(self):
		self.t.set_size(self.size)
		draw()

	def _snap_changed(self):
		self.t.set_snap(self.snap)
		draw()

	def _stretch_changed(self):
		self.t.set_stretch(self.stretch)
		draw()

	def _style_changed(self):
		self.t.set_style(self.style)
		draw()

	def _text_changed(self):
		self.t.set_text(str(self.text))
		draw()

	def _url_changed(self):
		self.t.set_url(self.url)
		draw()

	def _variant_changed(self):
		self.t.set_variant(self.variant)
		draw()

	def _verticalalignment_changed(self):
		self.t.set_verticalalignment(self.verticalalignment)
		draw()

	def _visible_changed(self):
		self.t.set_visible(self.visible)
		draw()

	def _weight_changed(self):
		self.t.set_weight(self.weight)
		draw()

	def _zorder_changed(self):
		self.t.set_zorder(self.zorder)
		draw()

	traits_view = etua.View(
		etua.Item("get_current_texts", show_label = False),
		etua.Item("t_list_max", style="readonly"),
		"t_list_index",
		etua.Group(
			etua.Group(
				'alpha',
				etua.Item('bbox_facecolor'),# style='custom'),
				'bbox_alpha',
				'clip_on',
				etua.Item('color'), #style='custom'),
				'horizontalalignment',
				'label',
				'name',
				'position_x',
				'position_y',
				'rotation',
				'size',
				'snap',
				'stretch',
				'style',
				'text',
				'variant',
				'verticalalignment',
				'visible',
				'weight',
				'zorder',
				label='text', dock='tab'),
			etua.Group(
				'animated',
				'clip_box',
				'clip_path',
				'family',
				'picker',
				'url',
				label='other', dock='tab'),
			layout='tabbed'
			),
		resizable = True,
		scrollable = True
		)

def settexts():
	SetupText().configure_traits()
