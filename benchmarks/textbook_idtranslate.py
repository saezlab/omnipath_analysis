#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Copyright 2019-2020 Saez Lab
#
# OmniPath2 analysis and figures suite
#
# Authors:
#
# Dénes Türei
# turei.denes@gmail.com
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  Website: https://omnipathdb.org/
#

import os
import collections
import itertools

import pypath.utils.mapping as mapping
import pypath.share.common as common


datadir = 'data'
infile = os.path.join(datadir, 'cytokine_receptor_from_book.tsv')
outfile = os.path.join(datadir, 'cytokine_receptor_book_translated.tsv')


nonstandard_names = {
    'IL-14': 'P40222',
}
missing_receptors = {
    'IL-18': 'Q13478',
    'IL-2': {'P01589', 'P14784', 'P31785'},
    'IL-20': {'Q6UXL0', 'Q9UHF4'},
    'IL-4': 'P24394',
    'CCL18': {'P51685', 'Q99527', 'Q9BZ71'},
    'CXCL4': {'P49682', 'P25024'},
}


class TranslateTextbookTable(object):


    def __init__(self, infile = None, outfile = None):

        self.infile = infile or globals()['infile']
        self.outfile = outfile or globals()['outfile']


    def main(self):

        self.load()
        self.translate()
        self.export()


    def load(self):

        with open(self.infile, 'r') as fp:

            columns = [
                col.lower().replace('.', '_')
                for col in fp.readline().strip().split('\t')
            ]
            self._raw_record = collections.namedtuple(
                'TextbookRawRecord',
                columns,
            )

            self.raw_data = [
                self._raw_record(*line.strip('\n\r').split('\t'))
                for line in fp
            ]


    def translate(self):

        self._record = collections.namedtuple(
            'TextbookRecord',
            self._raw_record._fields + (
                'ligand_uniprot',
                'receptor_uniprot',
                'ligand_genesymbol',
                'receptor_genesymbol',
                'record_id',
            ),
        )

        self.result = [
            self._record(*record, *translated, str(i))
            for i, record in enumerate(self.raw_data)
            for translated in self.translate_record(record)
        ]


    @classmethod
    def translate_record(cls, record):

        empty = {''}

        ligand_uniprots = cls.translate_one(
            record.source_ensembl,
            record.name,
            record.synonyms,
        ) or empty
        receptor_uniprots = (
            cls.translate_one(
                record.receptor_ensembl,
                record.cytokine_receptor,
            ) or
            (
                common.to_set(missing_receptors[record.name])
                    if record.name in missing_receptors else
                empty
            )
        )

        for u_lig, u_rec in itertools.product(
            ligand_uniprots,
            receptor_uniprots,
        ):

            g_lig = mapping.map_name0(u_lig, 'uniprot', 'genesymbol') or ''
            g_rec = mapping.map_name0(u_rec, 'uniprot', 'genesymbol') or ''

            yield u_lig, u_rec, g_lig, g_rec


    @classmethod
    def translate_one(cls, ensembl, genesymbol, synonyms = None):

        uniprots = mapping.map_name(ensembl, 'ensembl', 'uniprot')

        if not uniprots and genesymbol in nonstandard_names:

            uniprots = common.to_set(nonstandard_names[genesymbol])

        if not uniprots:

            uniprots = cls.translate_genesymbol(genesymbol)

        if not uniprots and synonyms:

            for syn in synonyms.split(','):

                uniprots = cls.translate_genesymbol(syn)

                if uniprots:

                    break

        return uniprots


    @staticmethod
    def translate_genesymbol(genesymbol):

        uniprots = mapping.map_name(genesymbol, 'genesymbol', 'uniprot')

        if not uniprots:

            uniprots = mapping.map_name(
                genesymbol.replace('-', ''),
                'genesymbol',
                'uniprot',
            )

        return uniprots


    def export(self):

        with open(self.outfile, 'w') as fp:

            _ = fp.write('\t'.join(self._record._fields))
            _ = fp.write(os.linesep)
            _ = fp.write(
                os.linesep.join(
                    '\t'.join(rec) for rec in self.result
                )
            )


if __name__ == '__main__':

    _ = TranslateTextbookTable(infile, outfile).main()