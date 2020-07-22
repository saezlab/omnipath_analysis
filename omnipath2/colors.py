#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Nicolàs Palacio-Escat
# nicolas.palacio@bioquant.uni-heidelberg.de
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: http://omnipathdb.org/
#

import os
import re
import importlib as imp
import collections

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
import colorsys

from pypath.share import session as session_mod

from omnipath2 import settings as op2_settings


recolor = re.compile(
    r'([0-9]+)\t'
    r'([0-9]+)\t'
    r'([0-9]+)\t'
    r'(\w+)'
)

class Colors(session_mod.Logger):


    def __init__(self):

        session_mod.Logger.__init__(self, name = 'op2.colors')

        self.palette_dir = op2_settings.get('palette_dir')
        self.palettes = {}
        self.read_palettes()


    def reload(self):
        """
        Reloads the module and updates the class instance.
        """

        modname = self.__class__.__module__
        mod = __import__(modname, fromlist = [modname.split('.')[0]])
        imp.reload(mod)
        new = getattr(mod, self.__class__.__name__)
        setattr(self, '__class__', new)


    def read_palettes(self):

        self._log('Reading palettes from directory `%s`.' % self.palette_dir)

        for fname in os.listdir(self.palette_dir):

            if os.path.splitext(fname)[-1] == '.gpl':

                self.read_palette(os.path.join(self.palette_dir, fname))


    def read_palette(self, path, name = None):

        self._log('Loading color palette from `%s`.' % path)

        colors = collections.OrderedDict()

        with open(path, 'r') as fp:

            for l in fp:

                if not name and l.startswith('Name:'):

                    name = l.split(':')[-1].strip()

                m = recolor.search(l)

                if not m:

                    continue

                r, g, b, label = m.groups(0)

                colors[label] = self.scale_to_01((r, g, b))

        if not name:

            name = os.path.basename(os.path.splitext(path)[0])

        self.palettes[name] = colors

        self._log(
            'Color palette `%s` with %u colors from `%s`.' % (
                name,
                len(colors),
                path,
            )
        )


    @staticmethod
    def scale_to_01(values):

        values = tuple(map(float, values))

        return (
            tuple(i / 255. for i in values)
                if any(v > 1. for v in values) else
            values
        )


    def get_palette(self, name):
        """
        Returns a palette as a ``matplotlib.colors.ListedColormap``
        object.
        """

        return mpl.colors.ListedColormap(
            list(self.palettes[name].values())
        )


    def get_altered_tone(self, name, idx, factor):
        """
        Retrieves a color `idx` from the palette `name` with a lightness
        altered by `factor`.
        """

        return self.lightness(self.palettes[name][idx], factor)


    def lightness_palette(self, name, factor):
        """
        name : str
            Palette name.
        """

        return [
            self.lightness(rgb)
            for rgb in self.palettes[name]
        ]


    @staticmethod
    def lightness(rgb, factor):
        """
        rgb : tuple
            Color as an RGB tuple.
        factor : float
            Lightness adjustment factor.
        """

        hls = colorsys.rgb_to_hls(*rgb)
        hls = (
            hls[0],
            min(hls[1] * factor, 1.),
            hls[2],
        )

        return colorsys.hls_to_rgb(*hls)


    def __repr__(self):

        return '<Colors: %u palettes loaded>' % len(self)


    def __len__(self):

        return len(self.palettes)


    def show(self, name):
        """
        Shows a palette (in case you want to see its colors.
        """

        pal = self.get_palette(name)
        self.show_palette(pal)


    @staticmethod
    def show_palette(palette):

        colors = palette.colors if hasattr(palette, 'colors') else palette

        fig, ax = plt.subplots(figsize = (1, 4), dpi = 100)

        for x, color in enumerate(colors):

            ax.scatter(x, 0, c = [color], s = 1000)

        ax.axes.set_aspect(8.0)
        fig.show()