"""IGVF MPRA standard file schemas.

Verbatim copies of the JSON Schemas published in **MPRAlib**
(Max Schubach, Berlin Institute of Health at Charité;
https://github.com/kircherlab/MPRAlib), MIT licence:

    Copyright (c) 2025, Max Schubach, Berlin Institute of Health at
    Charite - Universitatsmedizin Berlin

    Permission is hereby granted, free of charge, to any person
    obtaining a copy of this software and associated documentation
    files (the "Software"), to deal in the Software without
    restriction, including without limitation the rights to use,
    copy, modify, merge, publish, distribute, sublicense, and/or
    sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following
    conditions: The above copyright notice and this permission
    notice shall be included in all copies or substantial portions
    of the Software.

These describe the community file formats agreed by the IGVF MPRA
focus group (Rosen et al. 2025, Supplementary Note S1). They are
reproduced rather than paraphrased on purpose: they define an
interchange standard, so a retyped near-copy would be worse than
useless. Serialised as Python dicts so the validator needs no
package-data plumbing at runtime.
"""

from __future__ import annotations

from typing import Any

SCHEMAS: "dict[str, dict[str, Any]]" = {
    "reporter_sequence_design": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA reporter sequence design',
            'type': 'object',
            'description': 'Description of the MPRA design.',
            'properties': {   'name': {   'type': 'string',
                                          'description': 'A unique-within-file identifier, one '
                                                         'unique string per designed sequence',
                                          'minLength': 1},
                              'sequence': {   'type': 'string',
                                              'description': 'DNA string of designed sequence, '
                                                             'consisting of A,C,G,T and no N',
                                              'minLength': 1,
                                              'pattern': '^[ATGCatgc]+$'},
                              'category': {   'type': 'string',
                                              'description': 'Category of designed sequence',
                                              'enum': [   'variant',
                                                          'element',
                                                          'synthetic',
                                                          'scrambled']},
                              'class': {   'type': 'string',
                                           'description': 'Class of designed sequence',
                                           'enum': [   'test',
                                                       'variant positive control',
                                                       'variant negative control',
                                                       'element active control',
                                                       'element inactive control']},
                              'source': {   'type': 'string',
                                            'description': 'Free-form description of the '
                                                           'origin of the sequence'},
                              'ref': {   'type': 'string',
                                         'description': 'reference sequence, e.g. GRCh38'},
                              'chr': {   'name': 'chr',
                                         'type': 'string',
                                         'description': 'Reference chromosome or contig name'},
                              'start': {   'type': 'integer',
                                           'description': '0-based position of the left-most '
                                                          'position of sequence with respect '
                                                          'to the reference chromosome',
                                           'minimum': 0},
                              'end': {   'type': 'integer',
                                         'description': '1-based position of the right-most '
                                                        'position of sequence with respect to '
                                                        'the reference chromosome',
                                         'minimum': 1},
                              'strand': {   'type': 'string',
                                            'description': 'strand of sequence in reference',
                                            'enum': ['+', '-']},
                              'variant_class': {   'description': 'the class of the '
                                                                  'variant(s), allowing for '
                                                                  'multiple variants to be '
                                                                  'tested in one sequence '
                                                                  '(haplotypes)',
                                                   'anyOf': [   {   'type': 'array',
                                                                    'minItems': 1,
                                                                    'items': {   'type': 'string',
                                                                                 'enum': [   'SNV',
                                                                                             'indel']}},
                                                                {   'type': 'string',
                                                                    'enum': ['NA']}]},
                              'variant_pos': {   'description': '0-based position of the start '
                                                                'of the normalized '
                                                                'representation of the '
                                                                'variant(s). integer within '
                                                                '[0, len(sequence) - 1]',
                                                 'anyOf': [   {   'type': 'array',
                                                                  'minItems': 1,
                                                                  'items': {   'type': 'integer',
                                                                               'minimum': 0}},
                                                              {   'type': 'string',
                                                                  'enum': ['NA']}]},
                              'SPDI': {   'description': '0-based, validated SPDI '
                                                         'representation of the variant(s), '
                                                         'e.g. NC_000001.11:25253603:G:A',
                                          'anyOf': [   {   'type': 'array',
                                                           'minItems': 1,
                                                           'items': {   'type': 'string',
                                                                        'pattern': '^[A-Za-z0-9_.]+:[0-9]+:[A-Za-z]*:[A-Za-z]*$'}},
                                                       {'type': 'string', 'enum': ['NA']}]},
                              'allele': {   'description': 'the allele of the variant(s) with '
                                                           'respect to the referenece '
                                                           'chromosome sequence',
                                            'anyOf': [   {   'type': 'array',
                                                             'minItems': 1,
                                                             'items': {   'type': 'string',
                                                                          'enum': [   'ref',
                                                                                      'alt']}},
                                                         {'type': 'string', 'enum': ['NA']}]},
                              'info': {   'type': 'string',
                                          'description': 'any additional comment or '
                                                         'information',
                                          'items': {'type': 'string'}}},
            'required': [   'name',
                            'sequence',
                            'category',
                            'class',
                            'variant_class',
                            'variant_pos',
                            'SPDI',
                            'allele'],
            'additionalProperties': False},
    "reporter_barcode_to_element_mapping": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'Report Barcode to Element Mapping',
            'type': 'object',
            'description': 'Creates the link between tested oligos to associated barcodes. Can '
                           'be pre-designed or learned by association sequencing.',
            'properties': {   'barcode': {   'type': 'string',
                                             'description': 'Barcode sequence. Allowed chars '
                                                            '[A,T,G,C]',
                                             'minLength': 1,
                                             'pattern': '^[ATGC]+$'},
                              'oligoName': {   'type': 'string',
                                               'description': 'Name of the oligo barcode is '
                                                              'assigned to.',
                                               'minLength': 1}},
            'required': ['barcode', 'oligoName'],
            'additionalProperties': False},
    "reporter_experiment_barcode": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Experiment Barcode',
            'type': 'object',
            'description': 'This format is needed to save the complete measurement on a '
                           'barcode level of an experiment.',
            'properties': {   'barcode': {   'type': 'string',
                                             'description': 'Barcode, allowed chars [A,T,G,C]',
                                             'minLength': 1,
                                             'pattern': '^[ATGC]+$'},
                              'oligo_name': {   'type': 'string',
                                                'description': 'Name of the oligo of the '
                                                               'design.',
                                                'minLength': 1}},
            'patternProperties': {   '^dna_count_': {   'anyOf': [   {'type': 'integer'},
                                                                     {   'type': 'string',
                                                                         'maxLength': 0}]},
                                     '^rna_count_': {   'anyOf': [   {'type': 'integer'},
                                                                     {   'type': 'string',
                                                                         'maxLength': 0}]}},
            'minProperties': 4,
            'required': ['barcode', 'oligo_name'],
            'additionalProperties': False},
    "reporter_experiment": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Experiment Barcode',
            'type': 'object',
            'description': 'This format is needed to save the complete measurement on a '
                           'barcode level of an experiment.',
            'properties': {   'replicate': {   'type': 'string',
                                               'description': 'Name of the replicate.',
                                               'minLength': 1},
                              'oligo_name': {   'type': 'string',
                                                'description': 'Name of the oligo of the '
                                                               'design.',
                                                'minLength': 1},
                              'dna_counts': {   'type': 'integer',
                                                'description': 'Number of raw DNA counts.',
                                                'minimum': 0},
                              'rna_counts': {   'type': 'integer',
                                                'description': 'Number of raw RNA counts.',
                                                'minimum': 0},
                              'dna_normalized': {   'type': 'number',
                                                    'description': 'Number of '
                                                                   'normalized/scaled DNA '
                                                                   'counts (CPM), 4 decimals.'},
                              'rna_normalized': {   'type': 'number',
                                                    'description': 'Number of '
                                                                   'normalized/scaled RNA '
                                                                   'counts (CPM), 4 decimals.'},
                              'log2FoldChange': {   'type': 'number',
                                                    'description': 'Fold change (normalized '
                                                                   'rna/dna ratio, in log2 '
                                                                   'space),  4 decimals.'},
                              'n_bc': {   'type': 'integer',
                                          'description': 'Number of observed barcodes for the '
                                                         'oligo.',
                                          'minimum': 0}},
            'required': [   'replicate',
                            'oligo_name',
                            'dna_counts',
                            'rna_counts',
                            'dna_normalized',
                            'rna_normalized',
                            'log2FoldChange',
                            'n_bc'],
            'additionalProperties': False},
    "reporter_element": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Element',
            'type': 'object',
            'description': 'This format stores the raw format statistical activity analysis '
                           'for elements. It is dependent on a background/negative set '
                           'distribution.',
            'properties': {   'oligo_name': {   'type': 'string',
                                                'description': 'Name of tested oligo.',
                                                'minLength': 1},
                              'log2FoldChange': {   'type': 'number',
                                                    'description': 'Fold change (normalized '
                                                                   'output/input ratio, in '
                                                                   'log2 space).'},
                              'inputCount': {   'type': 'number',
                                                'description': 'Input count (DNA), normalized '
                                                               '(CPM), mean across '
                                                               'replicates.'},
                              'outputCount': {   'type': 'number',
                                                 'description': 'Output count (RNA), '
                                                                'normalized (CPM), mean across '
                                                                'replicates.'},
                              'minusLog10PValue': {   'type': 'number',
                                                      'description': '-log10 of P-value'},
                              'minusLog10QValue': {   'type': 'number',
                                                      'description': '-log10 of Q-value '
                                                                     '(FDR)'}},
            'required': [   'oligo_name',
                            'log2FoldChange',
                            'inputCount',
                            'outputCount',
                            'minusLog10PValue',
                            'minusLog10QValue'],
            'additionalProperties': False},
    "reporter_variant": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Variant',
            'type': 'object',
            'description': 'This format stores the raw format statistical activity analysis '
                           'for variants.',
            'properties': {   'variant_id': {   'type': 'string',
                                                'description': 'Variant ID in Canonical SPDI '
                                                               'format.',
                                                'pattern': '^[A-Za-z0-9_.]+:[0-9]+:[A-Za-z]*:[A-Za-z]*$'},
                              'log2FoldChange': {   'type': 'number',
                                                    'description': 'Fold change (alt '
                                                                   'output/input ratio divided '
                                                                   'by ref output/input ratio, '
                                                                   'in log2 space).'},
                              'inputCountRef': {   'type': 'number',
                                                   'description': 'Input count reference '
                                                                  'allele, normalized (CPM), '
                                                                  'mean across replicates.'},
                              'outputCountRef': {   'type': 'number',
                                                    'description': 'Output count reference '
                                                                   'allele, normalized (CPM), '
                                                                   'mean across replicates.'},
                              'inputCountAlt': {   'type': 'number',
                                                   'description': 'Input count alternative '
                                                                  'allele, normalized (CPM), '
                                                                  'mean across replicates.'},
                              'outputCountAlt': {   'type': 'number',
                                                    'description': 'Output count alternative '
                                                                   'allele, normalized (CPM), '
                                                                   'mean across replicates.'},
                              'minusLog10PValue': {   'type': 'number',
                                                      'description': '-log10 of P-value'},
                              'minusLog10QValue': {   'type': 'number',
                                                      'description': '-log10 of Q-value (FDR)'},
                              'postProbEffect': {   'type': 'number',
                                                    'description': 'Posterior probability of a '
                                                                   'regulatory effect.'},
                              'CI_lower_95': {   'type': 'number',
                                                 'description': 'Lower bound of a 95% interval '
                                                                'for the variant effect.'},
                              'CI_upper_95': {   'type': 'number',
                                                 'description': 'Upper bound of a 95% interval '
                                                                'for the variant effect'},
                              'variantPos': {   'type': 'integer',
                                                'description': '0-based position of the start '
                                                               'of the variant in the tested '
                                                               'sequence  -1 if aggregation of '
                                                               'multiple positions withins '
                                                               'tested sequences.',
                                                'minimum': -1},
                              'refAllele': {   'type': 'string',
                                               'description': 'Normalized Canonical SPDI '
                                                              'reference variant sequence, '
                                                              'allowed chars [A,T,G,C]. If '
                                                              'empty use 0.',
                                               'pattern': '^([ATGC]+|0)$'},
                              'altAllele': {   'type': 'string',
                                               'description': 'Normalized Canonical SPDI '
                                                              'alternative variant sequence, '
                                                              'allowed chars [A,T,G,C]. If '
                                                              'empty use 0.',
                                               'pattern': '^([ATGC]+|0)$'}},
            'required': [   'variant_id',
                            'log2FoldChange',
                            'inputCountRef',
                            'outputCountRef',
                            'inputCountAlt',
                            'outputCountAlt',
                            'minusLog10PValue',
                            'minusLog10QValue',
                            'postProbEffect',
                            'CI_lower_95',
                            'CI_upper_95',
                            'variantPos',
                            'refAllele',
                            'altAllele'],
            'additionalProperties': False},
    "reporter_genomic_element": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Genomic Element',
            'type': 'object',
            'description': 'Defines the activity of an element/region within a genome. Can '
                           'only be used when exact chromosome start and end location within a '
                           'reference genome is available. Cannot be used for shuffled or '
                           'modified elements.',
            'properties': {   'chrom': {   'type': 'string',
                                           'description': 'Reference sequence chromosome or '
                                                          'scaffold.'},
                              'chromStart': {   'type': 'integer',
                                                'description': 'Start position in chromosome, '
                                                               '0-based inclusive.',
                                                'minimum': 0},
                              'chromEnd': {   'type': 'integer',
                                              'description': 'End position in chromosome, '
                                                             '0-based exclusive.',
                                              'minimum': 0},
                              'name': {   'type': 'string',
                                          'description': 'Name of tested element or region.'},
                              'score': {   'type': 'number',
                                           'description': 'Indicates how dark the peak will be '
                                                          'displayed in the browser (0-1000).'},
                              'strand': {   'type': 'string',
                                            'description': '+ or - for strand, . for unknown',
                                            'enum': ['+', '-', '.']},
                              'log2FoldChange': {   'type': 'number',
                                                    'description': 'Fold change (normalized '
                                                                   'output/input ratio, in '
                                                                   'log2 space).'},
                              'inputCount': {   'type': 'number',
                                                'description': 'Input count (DNA), normalized '
                                                               '(CPM), mean across '
                                                               'replicates.'},
                              'outputCount': {   'type': 'number',
                                                 'description': 'Output count (RNA), '
                                                                'normalized (CPM), mean across '
                                                                'replicates.'},
                              'minusLog10PValue': {   'type': 'number',
                                                      'description': '-log10 of P-value'},
                              'minusLog10QValue': {   'type': 'number',
                                                      'description': '-log10 of Q-value '
                                                                     '(FDR)'}},
            'required': [   'chrom',
                            'chromStart',
                            'chromEnd',
                            'name',
                            'score',
                            'strand',
                            'log2FoldChange',
                            'inputCount',
                            'outputCount',
                            'minusLog10PValue',
                            'minusLog10QValue'],
            'additionalProperties': False},
    "reporter_genomic_variant": {   '$schema': 'http://json-schema.org/draft-07/schema#',
            'title': 'MPRA Reporter Genomic Variant',
            'type': 'object',
            'description': 'Defines the activity of a variant within a genome. Can only be '
                           'used when exact chromosome start and end location, reference and '
                           'alternative sequence (Canonical SPDI normalized format) within a '
                           'reference sequence is available. Cannot be used when no Canonical '
                           'SPDI is available for the variant.',
            'properties': {   'chrom': {   'type': 'string',
                                           'description': 'Reference sequence chromosome or '
                                                          'scaffold.'},
                              'chromStart': {   'type': 'integer',
                                                'description': 'Start position in chromosome, '
                                                               '0-based inclusive.',
                                                'minimum': 0},
                              'chromEnd': {   'type': 'integer',
                                              'description': 'End position in chromosome, '
                                                             '0-based exclusive.',
                                              'minimum': 0},
                              'name': {   'type': 'string',
                                          'description': 'Name of tested variant.'},
                              'score': {   'type': 'number',
                                           'description': 'Indicates how dark the peak will be '
                                                          'displayed in the browser (0-1000).'},
                              'strand': {   'type': 'string',
                                            'description': '+ or - for strand, . for unknown',
                                            'enum': ['+', '-', '.']},
                              'log2FoldChange': {   'type': 'number',
                                                    'description': 'Fold change (alt '
                                                                   'output/input ratio divided '
                                                                   'by ref output/input ratio, '
                                                                   'in log2 space).'},
                              'inputCountRef': {   'type': 'number',
                                                   'description': 'Input count reference '
                                                                  'allele, normalized (CPM), '
                                                                  'mean across replicates.'},
                              'outputCountRef': {   'type': 'number',
                                                    'description': 'Output count reference '
                                                                   'allele, normalized (CPM), '
                                                                   'mean across replicates.'},
                              'inputCountAlt': {   'type': 'number',
                                                   'description': 'Input count alternative '
                                                                  'allele, normalized (CPM), '
                                                                  'mean across replicates.'},
                              'outputCountAlt': {   'type': 'number',
                                                    'description': 'Output count alternative '
                                                                   'allele, normalized (CPM), '
                                                                   'mean across replicates.'},
                              'minusLog10PValue': {   'type': 'number',
                                                      'description': '-log10 of P-value'},
                              'minusLog10QValue': {   'type': 'number',
                                                      'description': '-log10 of Q-value (FDR)'},
                              'postProbEffect': {   'type': 'number',
                                                    'description': 'Posterior probability of a '
                                                                   'regulatory effect.'},
                              'CI_lower_95': {   'type': 'number',
                                                 'description': 'Lower bound of a 95% interval '
                                                                'for the variant effect.'},
                              'CI_upper_95': {   'type': 'number',
                                                 'description': 'Upper bound of a 95% interval '
                                                                'for the variant effect'},
                              'variantPos': {   'type': 'integer',
                                                'description': '0-based position of the start '
                                                               'of the variant in the tested '
                                                               'sequence  -1 if aggregation of '
                                                               'multiple positions withins '
                                                               'tested sequences.',
                                                'minimum': -1},
                              'refAllele': {   'type': 'string',
                                               'description': 'Normalized Canonical SPDI '
                                                              'reference variant sequence, '
                                                              'allowed chars [A,T,G,C]. If '
                                                              'empty use 0.',
                                               'pattern': '^([ATGC]+|0)$'},
                              'altAllele': {   'type': 'string',
                                               'description': 'Normalized Canonical SPDI '
                                                              'alternative variant sequence, '
                                                              'allowed chars [A,T,G,C]. If '
                                                              'empty use 0',
                                               'pattern': '^([ATGC]+|0)$'}},
            'required': [   'chrom',
                            'chromStart',
                            'chromEnd',
                            'name',
                            'score',
                            'strand',
                            'log2FoldChange',
                            'inputCountRef',
                            'outputCountRef',
                            'inputCountAlt',
                            'outputCountAlt',
                            'minusLog10PValue',
                            'minusLog10QValue',
                            'postProbEffect',
                            'CI_lower_95',
                            'CI_upper_95',
                            'variantPos',
                            'refAllele',
                            'altAllele'],
            'additionalProperties': False},
}

SCHEMA_NAMES = tuple(SCHEMAS)

