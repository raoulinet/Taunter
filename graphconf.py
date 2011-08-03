import enthought.traits.api as eta
import enthought.traits.ui.api as etua
import enthought.traits.ui.menu as etum
import enthought.traits.ui.wx.tree_editor as etuwt
import enthought.pyface.api as epa
import simplejson as sj
import os
import math


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
	"""
	cf
	"""

	return gcf()



def cAs():
	"""
	ca
	"""

	return getp(cf(), 'axes')



def cA(num = -1):
	"""
	cA
	"""

	return getp(cf(), 'axes')[num]



def ci():
	"""
	ci
	"""

	return gci()



def cl(line_number = -1):
	"""
	cl
	"""

	return gca().lines[line_number]



def cls(lower = 0, higher = None):
	"""
	cls
	"""

	if higher == None:
		higher = len(gca().lines)

	return gca().lines[lower : higher]



def get_info():
	"""
	get_info
	"""

	n = 0

	for i in gca().lines:
		print("# " + str(n)
		+ "\tmarker: " + str(i.get_marker())
		+ "\tline: " + str(i.get_linestyle())
		+ "\tcolor: " + str(i.get_color())
		+ "\tlen: " + str(len(i.get_xdata())))
		n += 1



def ct(text_number = -1):
	"""
	ct
	"""

	return gca().texts[text_number]



def ct(text_number = -1):
	"""
	ct
	"""

	return gca().texts



def cx():
	"""
	cx
	"""
	
	return cA().xaxis



def cy():
	"""
	cy
	"""
	
	return cA().yaxis



def cT():
	"""
	cT
	"""
	
	return ca().title



def cPs():
	"""
	"""

	return filter(lambda i: type(i) == Polygon, getp(ca(), 'children'))



def cP(num = -1):
	"""
	"""

	return filter(lambda i: type(i) == Polygon, getp(ca(), 'children'))[num]



def cEs():
	"""
	"""

	return filter(lambda i: type(i) == Ellipse, getp(ca(), 'children'))



def cE(num = -1):
	"""
	"""

	return filter(lambda i: type(i) == Ellipse, getp(ca(), 'children'))[num]



def colorize(palette = "fancy", offset = 0, period = None):
	"""
	to colorize lines depdening of layers
	"""

	if type(palette) == str:
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
	elif type(palette) == list:
		colors = palette
	else:
		colors = default_colors

	if period == None:
		period = len(colors)

	n = len(gca().lines)

	for i in range(n):
		tmp_col = colors[i%period + offset%(len(colors) - 1)]
		setp(cl(i), 'color', tmp_col, 'mfc', tmp_col, 'mec', tmp_col)

	draw()



def get_data(layer = -1):
	"""
	get_data
	"""

	return gca().lines[layer].get_data()



def rotate_data(theta = 0):
	"""
	rotate_data
	"""

	for i in gca().lines:
		_x1, _y1 = i.get_data()
		i.set_xdata(cos(theta) * _x1 + sin(theta) * _y1)
		i.set_ydata(-sin(theta) * _x1 + cos(theta) * _y1)

	draw()
	


def translate_data(displacement = 0):
	"""
	translate_data
	"""

	for i in gca().lines:
		_x1, _y1 = i.get_data()
		i.set_xdata(_x1 + displacement[0])
		i.set_ydata(_y1 + displacement[1])

	draw()



def rotate_data_by_mouse():
	"""
	rotate_data_by_mouse
	"""

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
	"""
	translate_data_by_mouse	
	"""

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
	"""
	get_index_list
	"""

	n = 0

	for i in gca().lines:
		print("# " + str(n)+ "\tmarker: " + str(i.get_marker()) + "\tline: " + str(i.get_linestyle()) + "\tcolor: " + str(i.get_color()) + "\tlen: " + str(len(i.get_xdata())))
		n += 1
	return len(gca().lines)



def slice_graph(num = None):

	source = gcf()

	x, y = source.axes[-1].lines[num].get_data()
	marker = source.axes[-1].lines[num].get_marker()
	markersize = source.axes[-1].lines[num].get_markersize()
	linestyle = source.axes[-1].lines[num].get_linestyle()
	linewidth = source.axes[-1].lines[num].get_linewidth()
	color = source.axes[-1].lines[num].get_color()
	label = source.axes[-1].lines[num].get_label()
	figure()
	plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)



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
		plot(x, y, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)



def copy_paste_graph(source_fig, drop_fig, xoffset = 0, xmul = 1, yoffset = 0, ymul = 1):

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
		plot(xmul*array(x) + xoffset, ymul*array(y) + yoffset, marker=marker, markersize=markersize, linestyle=linestyle, linewidth=linewidth, color=color, label=label)



def get_norm_angle():

    a = ginput()[0]

    return {
        "norm": sqrt(a[0]**2 + a[1]**2),
        "angle_rad": arctan2(a[1], a[0]),
        "angle_deg": arctan2(a[1], a[0])*180./pi
    }



def remove_first(layer = 0):

	gca().lines.pop(layer)
	draw()



def remove_last(layer = None):

	if layer == None:
		gca().lines.pop()
	else:
		gca().lines.pop(layer)
	draw()



def line_matrix(xstep = 1, ystep = 1, bottom = None, top = None):

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

    for j in intx:
	for k in inty:
		cc[k][j] = 0

    if len(gca().lines) == 1:
        offset = 1
    else:
        offset = 0

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
                    cc[round(j)][round(k)] = n + offset
                except:
                    print("oups, stopped @ " + str(i) + "/" + str(len(data)))

    return xx, yy, cc



################################################################################
#                                                                              #
# Different class of usefull graph.                                            #
# * graph(...)                                                                 #
# * (im | p)graph(...)                                                         #
#                                                                              #
################################################################################

class graph:
	"""
	Overlaod of figure() + plot() function
	a 2 in 1
	graph(...) to open and plot
	plot(...) to append a line in the graph
	"""
	
	def __init__(self, *args):
		figure()
		plot(*args)



class imgraph:
	"""
	Overlaod of figure() + imshow() function
	a 2 in 1
	imgraph(...) to open and plot
	imshow(...) to append a line in the graph
	"""

	def __init__(self, *args):
		figure()
		imshow(*args)



class pgraph:
	"""
	Overlaod of figure() + pcolor() function
	a 2 in 1
	colorgraph(...) to open and plot
	pcolormesh(...) to append a line in the graph
	"""

	def __init__(self, *args, **kwargs):
		figure()
		pcolormesh(*args, **kwargs)

class polargraph:
	"""
	Overlaod of figure() + pcolor() function
	a 2 in 1
	colorgraph(...) to open and plot
	pcolormesh(...) to append a line in the graph
	"""

	def __init__(self, *args, **kwargs):
		figure()
		polar()
		pcolormesh(*args, **kwargs)



################################################################################
#                                                                              #
# Importating and saving graph and data.                                       #
# * (save | open)JSONgraph(...)                                                #
# * openNanoQt(graph | polar_graph)(...)                                       #
#                                                                              #
################################################################################

def saveJSONgraph(file_name = None):
	"""
	Save the graph in JSON format
	make a bundle of all graph's properties.
	"""

	if file_name == None:
		print("Filename needed")
		return

	# Here is the bundle we plan to save
	bundle = {}

	# First, save figure parameters coming from gcf()
	bundle['figure'] = {} # place to save figure's properties

	bundle['figure']['alpha']       = getp(gcf(), 'alpha')
	bundle['figure']['animated']    = getp(gcf(), 'animated')
	bundle['figure']['edgecolor']   = getp(gcf(), 'edgecolor')
	bundle['figure']['facecolor']   = getp(gcf(), 'facecolor')
	# bundle['figure']['figheight']   = getp(gcf(), 'figheight')
	# bundle['figure']['figwidth']    = getp(gcf(), 'figwidth')
	bundle['figure']['frameon']     = getp(gcf(), 'frameon')
	bundle['figure']['label']       = getp(gcf(), 'label')
	bundle['figure']['visible']     = getp(gcf(), 'visible')
	bundle['figure']['zorder']      = getp(gcf(), 'zorder')

	# Now, save axes.
	# Axes composed with properties with lines (data) and texts (lables).
	# Need a loop to scan all of the stuffs in all of the axes.
	bundle['axes'] = [] # list of axes

	for current_axis in gcf().axes:

		# First of all, pack the general properties of the axes.
		bundle['axes'].append({})

		bundle['axes'][-1]['adjustable']     = getp(current_axis, 'adjustable')
		bundle['axes'][-1]['alpha']          = getp(current_axis, 'alpha')
		bundle['axes'][-1]['anchor']         = getp(current_axis, 'anchor')
		bundle['axes'][-1]['animated']       = getp(current_axis, 'animated')
		bundle['axes'][-1]['aspect']         = getp(current_axis, 'aspect')
		bundle['axes'][-1]['autoscale_on']   = getp(current_axis, 'autoscale_on')
		bundle['axes'][-1]['axis_bgcolor']   = getp(current_axis, 'axis_bgcolor')
		bundle['axes'][-1]['axisbelow']      = getp(current_axis, 'axisbelow')
		bundle['axes'][-1]['frame_on']       = getp(current_axis, 'frame_on')
		bundle['axes'][-1]['label']          = getp(current_axis, 'label')
		bundle['axes'][-1]['legend']         = getp(current_axis, 'legend')
		bundle['axes'][-1]['position']       = [current_axis._position._get_x0(), current_axis._position._get_y0(), current_axis._position._get_width(), current_axis._position._get_height()]
		bundle['axes'][-1]['title']          = getp(current_axis, 'title')
		bundle['axes'][-1]['visible']        = getp(current_axis, 'visible')
		bundle['axes'][-1]['xbound']         = getp(current_axis, 'xbound')
		bundle['axes'][-1]['xlabel']         = getp(current_axis, 'xlabel')
		bundle['axes'][-1]['xscale']         = getp(current_axis, 'xscale')
		bundle['axes'][-1]['ybound']         = getp(current_axis, 'ybound')
		bundle['axes'][-1]['ylabel']         = getp(current_axis, 'ylabel')
		bundle['axes'][-1]['yscale']         = getp(current_axis, 'yscale')
		bundle['axes'][-1]['zorder']         = getp(current_axis, 'zorder')

		# Now, pack the line.
		bundle['axes'][-1]['lines'] = []

		for current_line in current_axis.lines:
			
			bundle['axes'][-1]['lines'].append({})

			bundle['axes'][-1]['lines'][-1]['alpha']           = getp(current_line, 'alpha')
			bundle['axes'][-1]['lines'][-1]['animated']        = getp(current_line, 'animated')
			bundle['axes'][-1]['lines'][-1]['antialiased']     = getp(current_line, 'antialiased')
			bundle['axes'][-1]['lines'][-1]['color']           = getp(current_line, 'color')
			bundle['axes'][-1]['lines'][-1]['dash_capstyle']   = getp(current_line, 'dash_capstyle')
			bundle['axes'][-1]['lines'][-1]['dash_joinstyle']  = getp(current_line, 'dash_joinstyle')
			bundle['axes'][-1]['lines'][-1]['drawstyle']       = getp(current_line, 'drawstyle')
			bundle['axes'][-1]['lines'][-1]['label']           = getp(current_line, 'label')
			bundle['axes'][-1]['lines'][-1]['linestyle']       = getp(current_line, 'linestyle')
			bundle['axes'][-1]['lines'][-1]['linewidth']       = getp(current_line, 'linewidth')
			bundle['axes'][-1]['lines'][-1]['marker']          = getp(current_line, 'marker')
			bundle['axes'][-1]['lines'][-1]['markeredgecolor'] = getp(current_line, 'markeredgecolor')
			bundle['axes'][-1]['lines'][-1]['markeredgewidth'] = getp(current_line, 'markeredgewidth')
			bundle['axes'][-1]['lines'][-1]['markerfacecolor'] = getp(current_line, 'markerfacecolor')
			bundle['axes'][-1]['lines'][-1]['markersize']      = getp(current_line, 'markersize')
			bundle['axes'][-1]['lines'][-1]['solid_capstyle']  = getp(current_line, 'solid_capstyle')
			bundle['axes'][-1]['lines'][-1]['solid_joinstyle'] = getp(current_line, 'solid_joinstyle')
			bundle['axes'][-1]['lines'][-1]['visible']         = getp(current_line, 'visible')
			bundle['axes'][-1]['lines'][-1]['xdata']           = getp(current_line, 'xdata').tolist()
			bundle['axes'][-1]['lines'][-1]['ydata']           = getp(current_line, 'ydata').tolist()
			bundle['axes'][-1]['lines'][-1]['zorder']          = getp(current_line, 'zorder')

		# Then, pack the texts
		bundle['axes'][-1]['texts'] = []

		for current_text in current_axis.texts:
			
			bundle['axes'][-1]['texts'].append({})

			bundle['axes'][-1]['texts'][-1]['alpha']               = getp(current_text, 'alpha')
			bundle['axes'][-1]['texts'][-1]['animated']            = getp(current_text, 'animated')
			bundle['axes'][-1]['texts'][-1]['color']               = getp(current_text, 'color')
			bundle['axes'][-1]['texts'][-1]['family']              = getp(current_text, 'family')
			bundle['axes'][-1]['texts'][-1]['horizontalalignment'] = getp(current_text, 'horizontalalignment')
			bundle['axes'][-1]['texts'][-1]['label']               = getp(current_text, 'label')
			bundle['axes'][-1]['texts'][-1]['position']            = getp(current_text, 'position')
			bundle['axes'][-1]['texts'][-1]['rotation']            = getp(current_text, 'rotation')
			bundle['axes'][-1]['texts'][-1]['size']                = getp(current_text, 'size')
			bundle['axes'][-1]['texts'][-1]['stretch']             = getp(current_text, 'stretch')
			bundle['axes'][-1]['texts'][-1]['style']               = getp(current_text, 'style')
			bundle['axes'][-1]['texts'][-1]['text']                = getp(current_text, 'text')
			bundle['axes'][-1]['texts'][-1]['variant']             = getp(current_text, 'variant')
			bundle['axes'][-1]['texts'][-1]['verticalalignment']   = getp(current_text, 'verticalalignment')
			bundle['axes'][-1]['texts'][-1]['visible']             = getp(current_text, 'visible')
			bundle['axes'][-1]['texts'][-1]['weight']              = getp(current_text, 'weight')
			bundle['axes'][-1]['texts'][-1]['zorder']              = getp(current_text, 'zorder')

	target_file = file(file_name, "w")
	sj.dump(bundle, target_file)
	target_file.close()




def openJSONgraph(bundle = None):
	"""
	Open the graph from JSON format
	make a bundle of all graph's properties.
	"""

	figure()

	# First, set figure parameters
	
	try:
		setp(gcf(), 'alpha'    , bundle['figure']['alpha']    , 
			    'animated' , bundle['figure']['animated'] , 
			    'edgecolor', bundle['figure']['edgecolor'], 
			    'facecolor', bundle['figure']['facecolor'], 
			    # 'figheight', bundle['figure']['figheight'], 
			    # 'figwidth' , bundle['figure']['figwidth'] , 
			    'frameon'  , bundle['figure']['frameon']  , 
			    'label'    , bundle['figure']['label']    , 
			    'visible'  , bundle['figure']['visible']  , 
			    'zorder'   , bundle['figure']['zorder']   ) 
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
				setp(gca().lines[-1], 'alpha'          , current_line['alpha']          ,
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

	figure()    
	# setp(cA(), 'xlabel', bundle['x_label'])
	# setp(cA(), 'ylabel', bundle['y_label'])

	if bundle.has_key("curves"):

		for n in range(len(bundle["curves"])):
			
			if bundle["curves"][n] != None:
				
				curvex, curvey = hsplit(array(bundle["curves"][n]["data"]), 2)
				plot(curvex, curvey)
				print("Plotted the curve #" + str(n))

				print("pen_style " + curve_line[bundle["curves"][n]["options"]["pen_style"]])
				if curve_line[bundle["curves"][n]["options"]["pen_style"]] == "Solid":
					setp(cl(), 'ls', '-')
				else:
					setp(cl(), 'ls', '')

				if curve_line[bundle["curves"][n]["options"]["pen_style"]] == "-":
					setp(cl(), 'ls', '-')
				else:
					setp(cl(), 'ls', '')

				setp(cl(), 'lw', 1)
				setp(cl(), 'marker', 's')
				setp(cl(), 'ms', 3)
		colorize()
	else:
		print("No data to process")



def vestige_openNanoQtgraph(bundle, reversed = False):

    figure()    
    # setp(cA(), 'title', bundle['title']) 
    # setp(cA(), 'xlabel', bundle['x_label'])
    # setp(cA(), 'ylabel', bundle['y_label'])
    
    # if bundle['logscale_x']:
    #    setp(cA(), 'xscale', 'log')
    # else:
    #    setp(cA(), 'xscale', 'linear')
    
    # if bundle['logscale_y']:
    #    setp(cA(), 'yscale', 'log')
    # else:
    #    setp(cA(), 'yscale', 'linear')
    
    if bundle.has_key("curves"):

	print("Load the curves")

        for n in range(len(bundle["curves"])):

	    print("Seek data in the curve layer #" + str(n))

            if bundle["curves"][n] != None:

                curvex, curvey = hsplit(array(bundle["curves"][n]["data"]), 2)
                plot(curvex, curvey)
		print("Plotted the curve #" + str(n))

		# setp(cl(), 'ls', curve_line[bundle["curves"][n]["options"]["pen_style"]])

                # if str(bundle["curves"][n]["options"]["pen_color"]).startswith("0x"):
                #     formatted_color = '#' + str(bundle["curves"][n]["options"]["pen_color"])[2:]
                # else:
                #    formatted_color = str(bundle["curves"][n]["options"]["pen_color"])

		# setp(gca().lines[-1], 'c', formatted_color)
                # setp(gca().lines[-1], 'mec', formatted_color)
                # setp(gca().lines[-1], 'mfc', formatted_color)
    
                # setp(gca().lines[-1], 'marker', curve_marker[bundle["curves"][n]["options"]["symbol"]])
                # setp(gca().lines[-1], 'markersize', bundle["curves"][n]["options"]["symbol_size"])
    
    if bundle.has_key("array"):

    	print("Load array")

        a = array(bundle["array"]["data"]) 
        x, y = meshgrid(linspace(bundle["x_min"], bundle["x_max"], len(a)), linspace(bundle["y_min"], bundle["y_max"], len(a[0])))
        if reversed:
            pcolormesh(x, y, a.max() - a)
        else:
            pcolormesh(x, y, a)

    setp(cl(), 'ls', "")
    setp(cl(), 'lw', 1)
    colorize()

    draw()



def openNanoQtgraph_polar(bundle, reversed = False, deg_rad = 1):
    
    figure()    
    polar()
    if bundle.has_key("array"):
        a = array(bundle["array"]["data"]) 
        x, y = meshgrid(linspace(bundle["x_min"], bundle["x_max"], len(a)), linspace(bundle["y_min"], bundle["y_max"], len(a[0])))
        if reversed:
            pcolormesh(deg_rad*y, x, a.max() - a)
        else:
            pcolormesh(deg_rad*y, x, a)

    draw()



def openGraph(fname, reversed = False):

	if fname[-7:-1] == ".pyplo":
	
		try:
			data = sj.load(file(fname, "r"))
		except:
			print("Warning: fail to open the JSON graph")
			return
		
		try:
			data.has_key('axes')
		except:
			print("does not contains axes...")
			return

		try:
			openJSONgraph(data)
		except:
			print("Warning: fail to process the JSON graph")
			return

	elif fname[-5:-1] == ".plo":

		try:
			data = sj.load(file(fname, "r"))
			try:
				openNanoQtgraph(data, reversed = reversed)
			except:
				print("Warning: ouch @ processing!")
				return
		except:
			print("Warning: ouch @ opening!")
			return
	else:
		print("type not known")
		return



def openPolar(fname, reversed = False, angle = 2 * math.pi):

	try:
		data = sj.load(file(fname, "r"))
		try:
			openNanoQtgraph_polar(data, reversed = reversed, deg_rad = eval(angle))
		except:
			print("Warning: ouch @ processing!")
	except:
		print("Warning: ouch @ opening!")



class UiImport(eta.HasTraits):

    data = eta.Dict()
    fname = eta.File()
    angle = eta.Str("2*pi")
    reversed = eta.Bool()
    open_graph = eta.Button("Open graph")
    # open_NanoQt_graph = eta.Button("Open NanoQt graph")
    open_polar = eta.Button("Open polar")
    reload = eta.Button("Reload")
   
 
    def __init__(self, fname = ""):
        self.fname = _ip.magic("pwd ")

    def _reload_fired(self):
        self.fname = _ip.magic("pwd ")

    def _open_graph_fired(self):

	if self.fname[-7:-1] == ".pyplo":
	
		try:
		    self.data = sj.load(file(self.fname, "r"))
		except:
		    print("Warning: fail to open the JSON graph")
		    return
		
		try:
		    self.data.has_key('axes')
		except:
		    print("does not contains axes...")
		    return

		try:
		    openJSONgraph(self.data)
		except:
		    print("Warning: fail to open the JSON graph")
		    return
	else:

		try:
		    self.data = sj.load(file(self.fname, "r"))
		    try:
			openNanoQtgraph(self.data, reversed = self.reversed)
		    except:
			print("Warning: fail to open the NanoQt graph")
			return
		except:
		    print("Warning: fail to open the NanoQt object")
		    return


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

    # open_nanoqt_graph = etum.Action(name = 'open NanoQt graph', action = '_open_NanoQt_graph_fired')

    open_polar = etum.Action(name = 'open polar', action = '_open_polar_fired')
    reload = etum.Action(name = 'reload', action = '_reload_fired')


    view = etua.View(
        # etua.Item("data", style="simple"),
        etua.Item('angle'),
        etua.Item('reversed'),
        # etua.Item('fname', editor=etua.FileEditor(filter = ['*.plot'], auto_set = True), style = "custom"),
        etua.Item('fname', editor=etua.FileEditor(auto_set = True), style = "custom"),
        resizable = True,
        scrollable = True,
	title= "UiImport",
	toolbar = etum.ToolBar(open_graph, open_polar, reload),
        height = 640,
        width = 800
    )



def uiimport(fname = ""):

    wxuiimport = UiImport(fname)
    wxuiimport.configure_traits()



class UiBrowser(eta.HasTraits):

    directory = eta.Directory()
    current_directory = eta.Str("")
    cd = eta.Button("Directory")
   
 
    def __init__(self):
        self.directory = _ip.magic("pwd ")
        self.current_directory = self.directory


    def _cd_fired(self):

	_ip.magic("cd " + self.directory)
        self.current_directory = _ip.magic("pwd ")

    change_directory = etum.Action(name = 'Directory', action = '_cd_fired')

    view = etua.View(
        etua.Item('current_directory'),
        etua.Item('directory', editor=etua.DirectoryEditor(), style = "custom"),
        resizable = True,
        scrollable = True,
	title= "Browser",
	toolbar = etum.ToolBar(change_directory),
        height = 640,
        width = 800
    )



def browser():

    wxbrowser = UiBrowser()
    wxbrowser.configure_traits()
