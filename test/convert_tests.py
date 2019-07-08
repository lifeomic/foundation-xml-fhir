from mock import patch
from unittest import TestCase
from src.convert import process
import os.path
import filecmp

results_payload_dict = {
    'FinalReport': {
        'PMI': {
            'LastName': 'doe',
            'FirstName': 'jane',
            'MRN': '1234',
            'Gender': 'Female',
            'DOB': '1970-01-01',
            'CollDate': '2000-01-01',
            'SubmittedDiagnosis': 'Cancer'
        },
        'Sample': {
            'TestType': 'Test 1'
        }
    },
    'variant-report': {
        'samples': {
            'sample': {
                '@name': 'sample1'
            }
        },
        'short-variants': {
            'short-variant': [
                {
                    '@gene': 'gene1',
                    '@cds-effect': '229C&gt;A',
                    '@functional-effect': 'missense',
                    '@allele-fraction': 0.488,
                    '@position': 'chr1:100',
                    '@depth': 200,
                    '@transcript': 'NM_001',
                    '@status': 'known',
                    '@protein-effect': 'R77S',
                    '@strand': '-'
                }
            ]
        },
        'copy-number-alterations': {
            'copy-number-alteration': [
                {
                    '@gene': 'CDK4',
                    '@position': 'chr12:58093932-58188144',
                    '@copy-number': '44',
                    '@status': 'known',
                    '@type': 'amplification',
                    '@number-of-exons': '5 of 5'
                }
            ]
        },
        'rearrangements': {
            'rearrangement': [
                {
                    '@status': 'known',
                    '@targeted-gene': 'CDK4',
                    '@type': 'truncation',
                    '@pos1': 'ch17:29557687-29887856',
                    '@pos2': 'ch6:66426718-66427149'
                }
            ]
        },
        'biomarkers': {
            'microsatellite-instability': {
                '@status': 'MSS'
            },
            'tumor-mutation-burden': {
                '@status': 'low',
                '@score': '0.73',
                '@unit': 'mutations-per-megabase'
            }
        }
    }
}


class Args:
    pass

class ConvertTest(TestCase):
    def setUp(self):
        self.args = Args()
        self.args.project_id = 'project1'
        self.args.subject_id = 'subject1'
        self.args.fasta = 'genome.fasta'
        self.args.genes = 'genes.ref'
        self.args.file_url = None
        self.args.vcf_out_file = None
        self.args.sequence_id = 'sequence_id'

    @patch("src.convert.parse_hgvs")
    def test_convert_with_subject(self, mock_parse_hgvs):
        self.maxDiff = None
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'

        fhir_resources = process(results_payload_dict, self.args)

        mock_parse_hgvs.assert_called_once_with('NM_001:c.229C>A', self.args.fasta, self.args.genes)

        specimen = fhir_resources[0]
        self.assertDictEqual(specimen, {
            'resourceType': 'Specimen',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': 'project1'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                     {
                    'system': 'http://lifeomic.com/fhir/variant-source',
                    'code': 'Foundation'
                }
                ]
            },
            'identifier': [
                {
                    'value': 'sample1'
                }
            ],
            'subject': {
                'reference': 'Patient/subject1'
            },
            'id': specimen['id']
        })

        sequence = fhir_resources[1]
        observation = fhir_resources[3]
        copy_number_observation = fhir_resources[4]
        rearrangement = fhir_resources[5]
        microsatellite = fhir_resources[6]
        tumor = fhir_resources[7]

        self.assertDictEqual(observation, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                           'valueCodeableConcept': {'coding': [{'code': '1100',
                                                                'display': 'gene1',
                                                                'system': 'http://www.genenames.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                           'valueCodeableConcept': {'coding': [{'code': '48004-6',
                                                                'display': 'NM_001:c.229C>A',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeType',
                           'valueCodeableConcept': {'coding': [{'code': 'LL380-7',
                                                                'display': 'missense',
                                                                'system': 'http://snomed.info/sct'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeName',
                           'valueCodeableConcept': {'coding': [{'code': '48005-3',
                                                                'display': 'p.R77S',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAllelicFrequency',
                           'valueCodeableConcept': {'coding': [{'code': '81258-6',
                                                                'display': 0.488,
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                           'valueCodeableConcept': {'coding': [{'code': '48002-0',
                                                                'display': 'somatic',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                           'valueCodeableConcept': {'coding': [{'code': '48001-2',
                                                                'display': '100',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAChromosome',
                           'valueCodeableConcept': {'coding': [{'code': '47999-8',
                                                                'display': 'chr1',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsTotalReadDepth',
                           'valueCodeableConcept': {'coding': [{'code': '82121-5',
                                                                'display': 200,
                                                                'system': 'http://loinc.org'}]}},
                         {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsVariantReadCount',
                           'valueCodeableConcept': {'coding': [{'code': '82121-5',
                                                                'display': '98',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsTranscriptID',
                           'valueCodeableConcept': {'coding': [{'code': '51958-7',
                                                                'display': 'NM_001',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}}
                        ],
            'id': observation['id'],
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'short'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'identifier': [{'system': 'https://lifeomic.com/observation/genetic',
                          'value': 'chr1:100:A:T'}],
            'status': 'final',
            'subject': {'reference': 'Patient/subject1'},
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': 'known',
                    'display': 'Foundation - Known'
                    }
                ]
            },
            'code': {
                'coding': [
                {
                    'system': 'http://loinc.org',
                    'code': '55233-1',
                    'display': 'Genetic analysis master panel-- This is the parent OBR for the panel holding all of the associated observations that can be reported with a molecular genetics analysis result.'
                }
                ]
            }
        })

        self.assertDictEqual(copy_number_observation, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                           'valueCodeableConcept': {'coding': [{'code': '1100',
                                                                'display': 'CDK4',
                                                                'system': 'http://www.genenames.org'}]}},
                         {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                           'valueCodeableConcept': {'coding': [{'code': '48004-6',
                                                                'display': 'Amplification: CN=44',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                           'valueCodeableConcept': {'coding': [{'code': '48002-0',
                                                                'display': 'somatic',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                           'valueCodeableConcept': {'coding': [{'code': '48001-2',
                                                                'display': '58093932-58188144',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAChromosome',
                           'valueCodeableConcept': {'coding': [{'code': '47999-8',
                                                                'display': 'chr12',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsCopyNumberEvent',
                           'valueCodeableConcept': {'coding': [{'code': 'SO:0001019',
                                                                'display': 'Amplification',
                                                                'system': 'http://www.sequenceontology.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeName',
                           'valueCodeableConcept': {'coding': [{'code': '48005-3',
                                                                'display': 'Exons 5 of 5',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-copyNumber',
                           'valueCodeableConcept': {'coding': [{'code': 'copyNumber',
                                                                'display': '44',
                                                                'system': 'http://lifeomic.com'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}}
                            ],
            'id': copy_number_observation['id'],
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'copy-number'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'status': 'final',
            'subject': {'reference': 'Patient/subject1'},
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': 'known',
                    'display': 'Foundation - Known'
                    }
                ]
            },
            'code': {
                'coding': [
                {
                    'system': 'http://loinc.org',
                    'code': '55233-1',
                    'display': 'Genetic analysis master panel-- This is the parent OBR for the panel holding all of the associated observations that can be reported with a molecular genetics analysis result.'
                }
                ]
            }
        })

        self.assertDictEqual(rearrangement, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                           'valueCodeableConcept': {'coding': [{'code': '1100',
                                                                'display': 'CDK4',
                                                                'system': 'http://www.genenames.org'}]}},
                         {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                           'valueCodeableConcept': {'coding': [{'code': '48004-6',
                                                                'display': 'CDK4 Truncation',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                           'valueCodeableConcept': {'coding': [{'code': '48002-0',
                                                                'display': 'somatic',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                           'valueCodeableConcept': {'coding': [{'code': '48001-2',
                                                                'display': 'ch17:29557687-29887856 ch6:66426718-66427149',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}}
                            ],
            'id': rearrangement['id'],
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'rearrangement'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'status': 'final',
            'subject': {'reference': 'Patient/subject1'},
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': 'known',
                    'display': 'Foundation - Known'
                    }
                ]
            },
            'code': {
                'coding': [
                {
                    'system': 'http://loinc.org',
                    'code': '55233-1',
                    'display': 'Genetic analysis master panel-- This is the parent OBR for the panel holding all of the associated observations that can be reported with a molecular genetics analysis result.'
                }
                ]
            }
        })

        self.assertDictEqual(microsatellite, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}}],
            'id': microsatellite['id'],
            'effectiveDateTime': '2000-01-01',
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'microsatellite-instability'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'status': 'final',
            'subject': {'reference': 'Patient/subject1'},
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://loinc.org',
                    'code': 'LA14122-8',
                    'display': 'Stable'
                    }
                ]
            },
            'code': {
                'coding': [
                    {
                        'system': 'http://loinc.org',
                        'code': '81695-9',
                        'display': 'Microsatellite instability [Interpretation] in Cancer specimen Qualitative.'
                    }
                ]
            },
        })

        self.assertDictEqual(tumor, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}}],
            'id': tumor['id'],
            'effectiveDateTime': '2000-01-01',
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'tumor-mutation-burden'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'status': 'final',
            'subject': {'reference': 'Patient/subject1'},
            'code': {
                'coding': [
                {
                    'system': 'http://lifeomic.com/fhir/biomarker',
                    'code': 'TMB',
                    'display': 'Tumor Mutation Burden'
                }
                ]
            },
            'component': [
                {
                    'code': {
                        'coding': [
                            {
                                'system': "http://lifeomic.com/fhir/biomarker",
                                'code': 'TMB Status',
                                'display': 'TMB Status'
                            }
                        ]
                    },
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': "http://lifeomic.com/fhir/biomarker",
                                'code': 'TMB-L',
                                'display': 'low'
                            }
                        ]
                    }
                },
                {
                    'code': {
                        'coding': [
                            {
                                'system': "http://lifeomic.com/fhir/biomarker",
                                'code': "TMB Score",
                                'display': "TMB Score"
                            }
                        ]
                    },
                    'valueQuantity': {
                        'value': 0.73,
                        'unit': 'mutations-per-megabase'
                    }
                }
            ]
        })

        self.assertDictEqual(sequence, {
            'resourceType': 'Sequence',
            'type': 'dna',
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                              {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    }, {
                    'system': 'http://lifeomic.com/fhir/variant-source',
                    'code': 'Foundation'
                }]},
            'patient': {
                'reference': 'Patient/subject1'
            },
            'specimen': {
                'display': 'sample1',
                'reference': 'Specimen/{}'.format(specimen['id'])
            },
            'referenceSeq': {
                'genomeBuild': 'GRCh37'
            },
            'id': sequence['id'],
            'variant': [
                {
                    'reference': 'Observation/{}'.format(observation['id'])
                },
                {
                    'reference': 'Observation/{}'.format(copy_number_observation['id'])
                },
                {
                    'reference': 'Observation/{}'.format(rearrangement['id'])
                },
                {
                    'reference': 'Observation/{}'.format(microsatellite['id'])
                },
                {
                    'reference': 'Observation/{}'.format(tumor['id'])
                }
            ]
        })

        report = fhir_resources[2]
        self.assertDictEqual(report, {
            'resourceType': 'DiagnosticReport',
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                             {
                'system': 'http://lifeomic.com/fhir/source',
                'code': 'LifeOmic Task Service'
            }, {
                    'system': 'http://lifeomic.com/fhir/report-source',
                    'code': 'Foundation'
                }]},
            'extension': [
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/DiagnosticReport-geneticsAssessedCondition',
                    'valueReference': {
                        'reference': 'Cancer'
                    }
                },
                {'url': 'http://lifeomic.com/fhir/StructureDefinition/sequence-id',
                 'valueString': 'sequence_id'
                }
            ],
            'status': 'final',
            'code': {
                'text': 'Test 1'
            },
            'effectiveDateTime': '2000-01-01',
            'subject': {
                'reference': 'Patient/subject1'
            },
            'specimen': [{
                'display': 'sample1',
                'reference': 'Specimen/{}'.format(specimen['id'])
            }],
            'result': [
                {
                    'reference': 'Observation/{}'.format(observation['id'])
                },
                {
                    'reference': 'Observation/{}'.format(copy_number_observation['id'])
                },
                {
                    'reference': 'Observation/{}'.format(rearrangement['id'])
                },
                {
                    'reference': 'Observation/{}'.format(microsatellite['id'])
                },
                {
                    'reference': 'Observation/{}'.format(tumor['id'])
                }
            ],
            'id': report['id']
        })

    @patch("src.convert.parse_hgvs")
    def test_convert_with_no_subject(self, mock_parse_hgvs):
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'
        self.args.subject_id = None

        fhir_resources = process(results_payload_dict, self.args)
        self.assertEquals(len(fhir_resources), 9)
        subject = fhir_resources[0]
        self.assertDictEqual(subject, {
            'resourceType': 'Patient',
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'},
                             {
                'system': 'http://lifeomic.com/fhir/source',
                'code': 'LifeOmic Task Service'
            }]},
            'name': [{
                'use': 'official',
                'family': 'doe',
                'given': ['jane']
            }],
            'identifier': [{
                'type': {
                    'coding': [{
                        'system': 'http://hl7.org/fhir/v2/0203',
                        'code': 'MR'
                    }]
                },
                'value': '1234'
            }],
            'gender': 'female',
            'birthDate': '1970-01-01',
            'id': subject['id']
        })

        specimen = fhir_resources[1]
        self.assertEquals(specimen['subject']['reference'],
                          'Patient/{}'.format(subject['id']))
        sequence = fhir_resources[2]
        self.assertEquals(sequence['patient']['reference'],
                          'Patient/{}'.format(subject['id']))
        report = fhir_resources[3]
        self.assertEquals(report['subject']['reference'],
                          'Patient/{}'.format(subject['id']))
        observation = fhir_resources[4]
        self.assertEquals(observation['subject']['reference'],
                          'Patient/{}'.format(subject['id']))


    @patch("src.convert.parse_hgvs")
    def test_convert_with_vcf(self, mock_parse_hgvs):
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'
        self.args.vcf_out_file = './tmp/subject.vcf'

        fhir_resources = process(results_payload_dict, self.args)
        # should just create report resource
        self.assertEquals(len(fhir_resources), 6)
        self.assertEquals(fhir_resources[0]['resourceType'], 'DiagnosticReport')
        self.assertEquals(fhir_resources[1]['resourceType'], 'Observation')
        self.assertEquals(fhir_resources[2]['resourceType'], 'Observation')
        self.assertTrue(os.path.isfile('./unsorted.vcf'))
        self.assertTrue(filecmp.cmp('./unsorted.vcf', './test/data/expected.vcf', shallow=False))


    @patch("src.convert.parse_hgvs")
    def test_convert_with_no_variants(self, mock_parse_hgvs):
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'
        self.args.vcf_out_file = './tmp/subject.vcf'

        no_variants_payload_dict = {
            'FinalReport': {
                'PMI': {
                    'LastName': 'doe',
                    'FirstName': 'jane',
                    'MRN': '1234',
                    'Gender': 'Female',
                    'DOB': '1970-01-01',
                    'CollDate': '2000-01-01',
                    'SubmittedDiagnosis': 'Cancer'
                },
                'Sample': {
                    'TestType': 'Test 1'
                }
            },
            'variant-report': {
                'samples': {
                    'sample': {
                        '@name': 'sample1'
                    }
                },
                'short-variants': None,
                'copy-number-alterations': None
            }
        }


        fhir_resources = process(no_variants_payload_dict, self.args)
        # should just create report resource
        self.assertEquals(len(fhir_resources), 1)
        self.assertEquals(fhir_resources[0]['resourceType'], 'DiagnosticReport')
        self.assertTrue(os.path.isfile('./unsorted.vcf'))
        self.assertTrue(filecmp.cmp('./unsorted.vcf', './test/data/expected_none.vcf', shallow=False))



    @patch("src.convert.parse_hgvs")
    def test_convert_with_one_variant(self, mock_parse_hgvs):
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'
        self.args.vcf_out_file = './tmp/subject.vcf'

        no_variants_payload_dict = {
            'FinalReport': {
                'PMI': {
                    'LastName': 'doe',
                    'FirstName': 'jane',
                    'MRN': '1234',
                    'Gender': 'Female',
                    'DOB': '1970-01-01',
                    'CollDate': '2000-01-01',
                    'SubmittedDiagnosis': 'Cancer'
                },
                'Sample': {
                    'TestType': 'Test 1'
                }
            },
            'variant-report': {
                'samples': {
                    'sample': {
                        '@name': 'sample1'
                    }
                },
                'short-variants': {
                    'short-variant':
                        {
                            '@gene': 'gene1',
                            '@cds-effect': '229C&gt;A',
                            '@functional-effect': 'missense',
                            '@allele-fraction': 0.488,
                            '@position': 'chr1:100',
                            '@depth': 200,
                            '@transcript': 'NM_001',
                            '@status': 'known',
                            '@protein-effect': 'R77S',
                            '@strand': '-'
                        }
                },
                'copy-number-alterations': {
                    'copy-number-alteration':
                        {
                            '@gene': 'CDK4',
                            '@position': 'chr12:58093932-58188144',
                            '@copy-number': '44',
                            '@status': 'known',
                            '@type': 'amplification',
                            '@number-of-exons': '5 of 5'
                        }
                }
            }
        }


        fhir_resources = process(no_variants_payload_dict, self.args)
        # should just create report resource
        self.assertEquals(len(fhir_resources), 3)
        self.assertEquals(fhir_resources[0]['resourceType'], 'DiagnosticReport')
        self.assertTrue(os.path.isfile('./unsorted.vcf'))
        self.assertTrue(filecmp.cmp('./unsorted.vcf', './test/data/expected.vcf', shallow=False))
