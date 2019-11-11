#!/usr/bin/env python
#-*- coding: utf-8 -*-

#
#  This file is part of the `omnipath2` Python module
#
#  Copyright
#  2019
#  Heidelberg University, Uniklinik RWTH Aachen
#
#  File author(s): Dénes Türei (turei.denes@gmail.com)
#                  Nicolàs Palacio
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://pypath.omnipathdb.org/
#

from future.utils import iteritems

import importlib as imp
import os
import time
import copy
import itertools

import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.backends.backend_pdf
import matplotlib.backends.backend_cairo
from matplotlib import ticker
from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.backends.backend_cairo
from matplotlib.backends.backend_cairo import FigureCanvasCairo

from pypath import session_mod
from pypath import common

import omnipath2
from omnipath2 import settings as op2_settings
from omnipath2 import colors
import omnipath2.path


def is_opentype_cff_font(filename):
    """
    This is necessary to fix a bug in matplotlib:
    https://github.com/matplotlib/matplotlib/pull/6714
    Returns True if the given font is a Postscript Compact Font Format
    Font embedded in an OpenType wrapper.  Used by the PostScript and
    PDF backends that can not subset these fonts.
    """
    
    if os.path.splitext(filename)[1].lower() == '.otf':
        
        result = _is_opentype_cff_font_cache.get(filename)
        
        if result is None:
            
            with open(filename, 'rb') as fd:
                
                tag = fd.read(4)
            
            result = (tag == b'OTTO')
            _is_opentype_cff_font_cache[filename] = result
        
        return result
    
    return False


mpl.font_manager.is_opentype_cff_font = is_opentype_cff_font


class PlotBase(omnipath2.path.PathBase):
    
    text_elements = (
        'axis_label',
        'ticklabel',
        'xticklabel',
        'title',
        'legend_title',
        'legend_label',
        'annotation',
    )
    timestamp = time.strftime('%Y%m%d')
    
    
    def __init__(
            self,
            fname = None,
            fname_param = (),
            fname_timestamp = True,
            timestamp_override = None,
            filetype = 'pdf',
            figures_dir = None,
            dir_timestamp = True,
            xlab = None,
            ylab = None,
            width = None,
            height = None,
            grid_rows = 1,
            grid_cols = 1,
            grid_hratios = None,
            grid_wratios = None,
            font_family = None,
            font_style = None,
            font_weight = None,
            font_variant = None,
            font_stretch = None,
            font_size = None,
            font_sizes = None,
            axis_label_font = None,
            ticklabel_font = None,
            xticklabel_font = None,
            legend_label_font = None,
            legend_title_font = None,
            title_font = None,
            annotation_font = None,
            legend_loc = 2,
            legend_title = None,
            legend = True,
            bar_args = None,
            xticks = None,
            xticklabels = None,
            xlim = None,
            ylim = None,
            uniform_ylim = False,
            palette = None,
            lab_size = (9, 9),
            axis_label_size = 10.0,
            xticklab_angle = 0,
            rc = None,
            title = None,
            maketitle = False,
            title_halign = None,
            usetex = False,
            do_plot = True,
            make_plot_first = False,
            log_label = None,
            plot_args = None,
            tight_layout = True,
            **kwargs,
        ):
        """
        Base class for various plotting classes. The ``plot`` method is the
        main method in this class, it executes the entire workflow of
        plotting. By default ``__init`` calls this method which means plot is
        made upon instantiation of the class. Plotting consists of three major
        phases done by methods ``pre_plot``, ``make_plot`` and ``post_plot``.
        In this base class ``pre_plot`` does nothing, ``make_plot`` creates
        a objects for figure, axes, grid, etc. It also does the plotting
        itself and applies the additional options like typeface, labels,
        sizes, etc. Then ``post_plot`` applies ``tight_layout`` and closes
        the file.
        This class is designed to be very much typography aware. You have
        easy control over the typefaces of all text elements of the figure.
        
        Parameters
        ----------
        fname : str
            File name for the graphics. Either a key for settings or a string
            provided directly e.g. `boxplot01`. The extension and optionally
            a timestamp will be added automatically. It might contain fields
            for formatting, the values for this can be provided in
            `fname_param`.
        fname_param : tuple
            Values to be inserted into the file name.
        fname_timestamp : bool
            Include timestamp in the file name. The timestamp will be the
            date of module load day in format of *YYYYMMDD*. It will be added
            before the file type extension, separated by double underscores.
            You can use different timestamp by providing it to the
            ``timestamp_override`` parameter.
        timestamp_override : str
            If provided overrides the default timestamp.
        filetype : str
            The graphics format to use. Default is *pdf*, alternatives are
            *png* or *svg*. Should correspond to the file name extension.
        figures_dir : str
            A path to the output directory for graphics files.
        xlab : str
            Label for the x axis.
        ylab : str
            Label for the y axis.
        width : float
            Figure width in inches.
        height : float
            Figure height in inches.
        figsizes : tuple
            Figure size as tuple of 2 numbers, in inches.
        grid_rows : int
            Number of rows in the grid. We create a grid even for single plot,
            hence the default is 1.
        grid_cols : int
            Number of columns in the grid. Just like ``grid_rows``.
        grid_hratios : list
            Height ratios of grid columns. List of floats.
        grid_wratios : list
            Width ratios of grid rows. List of floats.
        font_family : str,list
            Font family to use or list of families in order of preference.
            Default is from ``op2_settings.font_family``.
        font_style : str
            Font style, e.g. `bold` or `medium`.
        font_variant : str
            Font variant, e.g. `small-caps`.
        font_stretch : str
            Font stretch, e.g. `condensed` or `expanded`.
        font_size : float
            Base font size.
        font_sizes : dict
            Font size proportions. E.g. ``{'axis_label': 0.8}`` will
            result the axis label to have size of <base size> * 0.8.
        axis_label_font : dict
            Dict with the same font parameters as listed above, specific
            for axis labels.
        ticklabel_font : dict
            Dict with the same font parameters as listed above, specific
            for tick labels.
        xticklabel_font : dict
            Dict with the same font parameters as listed above, specific
            for tick labels of x axis.
        legend_title_font : dict
            Dict with the same font parameters as listed above, specific
            for the legend title.
        legend_label_font : dict
            Dict with the same font parameters as listed above, specific
            for the legend labels.
        title_font : dict
            Dict with the same font parameters as listed above, specific
            for the main title.
        annotation_font : dict
            Dict with the same font parameters as listed above, specific
            for the annotations.
        legend_loc : int
            The location of the legend.
        legend_title : str
            Title for the legend.
        legend : bool
            Create legend for the plot.
        bar_args : dict
            Arguments for barplots ``bar`` method.
        xticks : list
            Locations of ticks on x axis. List of floats.
        xticklabels : list
            Tick labels on x axis. List of strings.
        uniform_ylim : bool
            In case of multiple plots on a grid, the y axis limits should be
            uniform or different for each plot.
        palette : list
            Colours to use.
        lab_size : tuple
            Font size of the tick labels. Tuple of two numbers, for x and y
            axes, respectively.
        axis_label_size : float
            Font size of the axis labels.
        xticklab_angle : float
            Angle of the x axis labels.
        rc : dict
            Matplotlib rc params.
        title : str
            Main title of the plot.
        maketitle : bool
            Plot with or without main title.
        title_halign : str
            Horizontal alignement of the main title.
        usetex : bool
            Use LaTeX for rendering text elements.
        do_plot : bool
            Execute the plotting workflow upon instatiation.
            This is convenient by default but for resolving issues sometimes
            beneficial to first create the object and call individual methods
            afterwards.
        """
        
        if not hasattr(self, '_logger'):
            
            session_mod.Logger.__init__(self, name = log_label or 'op2.plot')
        
        figures_dir = figures_dir or omnipath2.data.figures_dir
        
        for k, v in itertools.chain(iteritems(locals()), iteritems(kwargs)):
            
            # we check this because derived classes might have set
            # already attributes
            if not hasattr(self, k) or getattr(self, k) is None:
                
                setattr(self, k, v)
        
        omnipath2.path.PathBase.__init__(
            self,
            fname = self.fname,
            fname_param = self.fname_param,
            fname_timestamp = self.fname_timestamp,
            timestamp_override = self.timestamp_override,
            filetype = self.filetype,
            target_dir = self.figures_dir,
        )
        
        if self.do_plot:
            
            self.main()
    
    
    def reload(self):
        """
        Reloads the module and updates the class instance.
        """
        
        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)
    
    
    def plot(self):
        """
        The total workflow of this class.
        Calls all methods in the correct order.
        """
        
        self.load_data()
        self.pre_plot()
        
        if self.make_plot_first:
            
            self.make_plot()
            self.set_figure()
            
            for ax in self.fig.axes:
                
                self.ax = ax
                self.post_subplot_hook()
            
        else:
            
            self.set_figure()
            self.make_plot()
        
        self.post_plot()
    
    # a synonym
    main = plot
    
    
    def load_data(self):
        
        pass
    
    
    def pre_plot(self):
        """
        Executes all necessary tasks before plotting in the correct order.
        Derived classes should override this if necessary.
        """
        
        omnipath2.path.PathBase.main(self)
        
        self._log('Plotting to `%s`.' % self.path)
        
        self.set_fonts()
        self.set_palette()
        self.set_rc()
        self.set_bar_args()
        self.set_plot_args()
        self.set_figsize()
    
    
    def set_figure(self):
        
        self.init_figure()
        self.set_grid()
    
    
    def make_plot(self):
        """
        Calls the plotting methods in the correct order.
        """
        
        self.make_plots() # this should call post_subplot_hook
                          # after each subplot
    
    
    def post_subplot_hook(self):
        """
        A method executed each time after creating a subplot.
        """
        
        self.set_labels()
        self.set_limits()
        self.set_ylims()
        self.set_title()
        self.set_ticklabels()
        self.make_legend()
        self.set_legend_font()
    
    
    def post_plot(self):
        """
        Saves the plot into file, and closes the figure.
        """
        
        self.finish()
    
    
    def set_bar_args(self):
        
        self.bar_args = self.bar_args or {}
    
    
    def set_plot_args(self):
        
        args = {
            'cmap': self.palette,
        }
        
        args.update(self.plot_args or {})
        
        self.plot_args = args
    
    
    def set_fonts(self):
        """
        Sets up everything related to fonts.
        """
        
        # if no font parameters provided by arguments or derived classes
        # we get the defaults from settings
        self.fonts_defaults_from_settings()
        # we create a dict of default font parameters
        self.fonts_default_dict()
        # replace None defaults with dicts
        self.fonts_init_dicts()
        # set font parameters in rc (is this necessary at all?)
        self.fonts_set_rc()
        # create font properties objects from all parameters dicts
        self.fonts_create_fontproperties()
    
    
    def fonts_init_dicts(self):
        """
        Initializes specific font argument dicts unless they are explicitely
        provided.
        """
        
        font_sizes = copy.deepcopy(op2_settings.get('font_sizes'))
        font_sizes.update(self.font_sizes or {})
        
        for text in self.text_elements:
            
            attr = '%s_font' % text
            this_font = copy.deepcopy(self.font_default)
            specific  = getattr(self, attr) or {}
            
            if 'size' not in specific:
                
                specific['size'] = self.font_size * (
                    font_sizes[text] if text in font_sizes else 1.0
                )
            
            this_font.update(specific)
            
            setattr(self, attr, this_font)
    
    
    def fonts_set_rc(self):
        """
        Sets up font related settings in matplotlib rc dict.
        """
        
        self.rc = self.rc or {}
        
        if type(self.lab_size) is not tuple:
            self.lab_size = (self.lab_size, ) * 2
        if 'axes.labelsize' not in self.rc:
            self.rc['axes.labelsize'] = self.axis_label_size
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[0]
        if 'ytick.labelsize' not in self.rc:
            self.rc['ytick.labelsize'] = self.lab_size[1]
        
        self.rc['font.family'] = self.font_family
        self.rc['font.style'] = self.font_style
        self.rc['font.variant'] = self.font_variant
        self.rc['font.weight'] = self.font_weight
        self.rc['font.stretch'] = self.font_stretch
        self.rc['text.usetex'] = self.usetex
    
    
    def fonts_defaults_from_settings(self):
        """
        Sets default font options from ``settings`` unless they are
        explicitely set already.
        """
        
        self.font_family = self.font_family or op2_settings.get('font_family')
        self.font_style = self.font_style or op2_settings.get('font_style')
        self.font_variant = (
            self.font_variant or op2_settings.get('font_variant')
        )
        self.font_stretch = (
            self.font_stretch or op2_settings.get('font_stretch')
        )
        self.font_size = self.font_size or op2_settings.get('font_size')
        self.font_weight = self.font_weight or op2_settings.get('font_weight')
    
    
    def fonts_default_dict(self):
        """
        Creates a dict from default font parameters which can be passed later
        as arguments for ``matplotlib.font_manager.FontProperties``.
        """
        
        # dict for default font parameters
        self.font_default = {
            'family': self.font_family,
            'style': self.font_style,
            'variant': self.font_variant,
            'stretch': self.font_stretch,
            'size': self.font_size,
            'weight': self.font_weight,
        }
    
    
    def fonts_create_fontproperties(self):
        """
        Creates ``matplotlib.font_manager.FontProperties`` objects from
        each font parameter dict.
        """
        
        self.fp = mpl.font_manager.FontProperties(
            **copy.deepcopy(self.font_default)
        )
        
        self.fp_default = (
            mpl.font_manager.FontProperties(
                family = self.font_family,
                style = self.font_style,
                weight = self.font_weight,
                variant = self.font_variant,
                stretch = self.font_stretch,
                size = self.font_size,
            )
        )
        
        
        for elem in self.text_elements:
            
            dictattr = '%s_font' % elem
            fpattr   = 'fp_%s' % elem
            
            fp = mpl.font_manager.FontProperties(
                **copy.deepcopy(getattr(self, dictattr))
            )
            
            setattr(self, fpattr, fp)
            
            self._log(
                'Font properties for `%s`: %s' % (
                    elem,
                    fp.get_fontconfig_pattern(),
                )
            )
            self._log(
                'Font file for `%s`: %s' % (
                    elem,
                    mpl.font_manager.fontManager.findfont(fp),
                )
            )
    
    def set_palette(self):
        
        self.palette = (
            self.palette
                if isinstance(self.palette, mpl.colors.Colormap) else
            colors.get_palette(self.palette)
                if self.palette in colors.palettes else
            mpl.cm.get_cmap(self.palette)
                if isinstance(self.palette, common.basestring) else
            mpl.colors.ListedColormap(self.palette)
                if isinstance(
                    self.palette,
                    (tuple, list, type({}.values()))
                ) else
            colors.get_palette(op2_settings.get('palette'))
        )
        
        mpl.cm.register_cmap(name = 'op2_current', cmap = self.palette)
        rcParams['image.cmap'] = 'op2_current'
    
    
    def set_rc(self):
        
        rcParams.update(self.rc or {})
    
    
    def set_figsize(self):
        """
        Converts width and height to a tuple so can be used for figsize.
        """
        
        if hasattr(self, 'figsize') and self.figsize is not None:
            
            return
            
        self.width = self.width or op2_settings.get('fig_width')
        self.height = self.height or op2_settings.get('fig_height')
        self.figsize = (self.width, self.height)
    
    
    def init_figure(self):
        """
        Creates a figure using the object oriented matplotlib interface.
        """
        
        self._create_figure()
        
        if self.filetype == 'pdf':
            
            self.pdf = mpl.backends.backend_pdf.PdfPages(self.path)
            self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        
        elif self.filetype in {'png', 'svg'}:
            
            self.cvs = (
                mpl.backends.backend_cairo.FigureCanvasCairo(self.fig)
            )
    
    
    def _create_figure(self):
        
        if not (
            hasattr(self, 'fig') and
            isinstance(self.fig, mpl.figure.Figure)
        ):
            
            self.fig = mpl.figure.Figure()
        
        self.fig.set_size_inches(self.figsize)
    
    
    def set_grid(self):
        """
        Sets up a grid according to the number of subplots,
        with proportions according to the number of elements
        in each subplot.
        """
        
        self.grid_hratios = self.grid_hratios or [1.] * self.grid_rows
        self.grid_wratios = self.grid_wratios or [1.] * self.grid_cols
        
        self.gs = mpl.gridspec.GridSpec(
            self.grid_rows,
            self.grid_cols,
            height_ratios = self.grid_hratios,
            width_ratios = self.grid_wratios,
        )
        
        self.axes = [[None] * self.grid_cols for _ in range(self.grid_rows)]
    
    
    def get_subplot(self, i = 0, j = 0, subplot_args = None):
        
        if self.axes[i][j] is None:
            
            subplot_args = subplot_args or {}
            
            self.axes[i][j] = self.fig.add_subplot(
                self.gs[i, j],
                **subplot_args
            )
        
        self.ax = self.axes[i][j]
    
    
    def add_subplot(self, ax, i = 0, j = 0):
        
        self.axes[i][j] = ax
    
    
    def iter_subplots(self):
        
        for j in xrange(self.grid_rows):
            
            for i in xrange(self.grid_cols):
                
                self.get_subplot(i, j)
                
                yield self.ax
    
    
    def make_plots(self):
        """
        By default this plots nothing here in the base class.
        Derived classes should override.
        """
        
        self.ax = get_subplot(0, 0)
        
        self.ax.plot(x = [], y = [])
        
        self.post_subplot_hook()
    
    
    def set_limits(self):
        
        if self.xlim:
            
            self.ax.set_xlim(self.xlim)
        
        if self.ylim:
            
            self.ax.set_ylim(self.ylim)
    
    
    def set_ylims(self):
        
        if self.uniform_ylim:
            
            maxy = max(
                ax.get_ylim()[1]
                for ax in self.axes[0][0:]
            )
            
            _ = [None for _ in ax.set_ylim([0, maxy]) for ax in self.axes[0]]
    
    
    def set_title(self):
        """
        Sets the main title.
        """
        
        if self.maketitle:
            
            self.title_text = self.fig.suptitle(self.title)
            self.title_text.set_fontproperties(self.fp_title)
            
            if self.title_halign:
                
                self.title_text.set_horizontalalignment(self.title_halign)
    
    
    def set_labels(self):
        """
        Sets properties of axis labels and ticklabels.
        """
        
        if self.xlab is not None:
            self._xlab = self.ax.set_xlabel(self.xlab)
        
        if self.ylab is not None:
            self._ylab = self.ax.set_ylabel(self.ylab)
        
        _ = [
            tick.label.set_fontproperties(self.fp_xticklabel) or (
                self.xticklab_angle == 0 or self.xticklab_angle == 90
            ) and (
                tick.label.set_rotation(self.xticklab_angle) or
                tick.label.set_horizontalalignment('center')
            )
            for tick in self.ax.xaxis.get_major_ticks()
        ]
        
        _ = [
            tick.label.set_fontproperties(self.fp_ticklabel)
            for tick in self.ax.yaxis.get_major_ticks()
        ]
        
        self.ax.set_ylabel(self.ylab, fontproperties = self.fp_axis_label)
        self.ax.set_xlabel(self.xlab, fontproperties = self.fp_axis_label)
        # self.ax.yaxis.label.set_fontproperties(self)
    
    
    def set_ticklabels(self):
        
        if self.xticklabels is not None:
            
            self.ax.set_xticklabels(
                ('{0}'.format(x) for x in self.xticklabels),
            )
        
        _ = [
            tl.set_fontproperties(self.fp_ticklabel)
            for tl in self.ax.get_xticklabels()
        ]
        
        if self.xticks is not None:
            
            self.ax.set_xticks(self.xticks)
    
    
    def make_legend(self):
        
        if self.legend:
            
            self.leg = self.ax.legend(
                loc = self.legend_loc,
                title = self.legend_title,
            )
    
    
    def set_legend_font(self):
        
        if self.legend:
            
            _ = [
                t.set_fontproperties(self.fp_legend_label)
                for t in self.leg.get_texts()
            ]
            
            self.leg.get_title().set_fontproperties(self.fp_legend_title)
    
    
    def finish(self):
        """
        Applies tight layout, draws the figure, writes the file and closes.
        """
        
        if self.tight_layout:
            
            self.fig.tight_layout()
        
        self.fig.subplots_adjust(
            top = (
                (.86 if '\n' in self.title else .9)
                    if self.maketitle else
                .92
            )
        )
        self.cvs.draw()
        
        self.save()
    
    
    def save(self):
        
        if self.filetype == 'pdf':
            
            self.cvs.print_figure(self.pdf)
            self.pdf.close()
            self._log('Plot saved to `%s`.' % self.path)
            
        elif self.filetype in {'png', 'svg'}:
            
            if self.path:
                
                with open(self.path, 'wb') as fp:
                    
                    getattr(self.cvs, 'print_%s' % self.filetype)(fp)
        
        self.ready()
        
        self._log('Plot saved to `%s`.' % self.path)
        
        #self.fig.clf()


class _TestPlot(PlotBase):
    """
    Most minimal class to test the ``PlotBase`` class with plotting only
    a sinus function. Also serves as an example how to build plot classes
    on top of PlotBase.
    """
    
    
    def __init__(self, **kwargs):
        
        defaults = {
            'fname': 'test_figure',
        }
        kwargs.update(defaults)
        
        PlotBase.__init__(self, **kwargs)
    
    
    def make_plots(self):
        
        self.get_subplot(0, 0)
        
        x = np.linspace(0, 10, 1000)
        self.ax.plot(x, np.sin(x))
        self.post_subplot_hook()
