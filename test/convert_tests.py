from mock import patch
from unittest import TestCase
from src.convert import process

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
                    '@protein-effect': 'R77S'
                }
            ]
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

        self.assertDictEqual(observation, {
            'extension': [{'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                           'valueCodeableConcept': {'coding': [{'code': '1100',
                                                                'display': 'gene1',
                                                                'system': 'http://www.genenames.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                           'valueCodeableConcept': {'coding': [{'code': '48004-6',
                                                                'display': '229C>A',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeType',
                           'valueCodeableConcept': {'coding': [{'code': 'LL380-7',
                                                                'display': 'missense',
                                                                'system': 'http://snomed.info/sct'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeName',
                           'valueCodeableConcept': {'coding': [{'code': '48005-3',
                                                                'display': 'R77S',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAllelicFrequency',
                           'valueCodeableConcept': {'coding': [{'code': '81258-6',
                                                                'display': 0.488,
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                           'valueCodeableConcept': {'coding': [{'code': '48002-0',
                                                                'display': 'somatic',
                                                                'system': 'http://loinc.org'}]}},
                          {'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                           'valueReference': {'reference': 'Sequence/{}'.format(sequence['id'])}},
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
                          {'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsTranscriptID',
                           'valueCodeableConcept': {'coding': [{'code': '51958-7',
                                                                'display': 'NM_001',
                                                                'system': 'http://loinc.org'}]}}],
            'id': observation['id'],
            'meta': {'tag': [{'code': 'project1',
                              'system': 'http://lifeomic.com/fhir/dataset'}]},
            'resourceType': 'Observation',
            'specimen': {'display': 'sample1',
                         'reference': 'Specimen/{}'.format(specimen['id'])},
            'identifier': [{'system': 'https://lifeomic.com/observation/genetic',
                          'value': 'chr1:100:A:T'}],
            'status': 'known',
            'subject': {'reference': 'Patient/subject1'}
        })

        self.assertDictEqual(sequence, {
            'resourceType': 'Sequence',
            'type': 'dna',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': 'project1'
                    }
                ]
            },
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
                }
            ]
        })

        report = fhir_resources[2]
        self.assertDictEqual(report, {
            'resourceType': 'DiagnosticReport',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': 'project1'
                    }
                ]
            },
            'extension': [
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/DiagnosticReport-geneticsAssessedCondition',
                    'valueReference': {
                        'reference': 'Cancer'
                    }
                }
            ],
            'status': 'final',
            'code': {
                'test': 'Test 1'
            },
            'issued': '2000-01-01',
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
                }
            ],
            'id': report['id']
        })

    @patch("src.convert.parse_hgvs")
    def test_convert_with_no_subject(self, mock_parse_hgvs):
        mock_parse_hgvs.return_value = 'chr1', 100, 'A', 'T'
        self.args.subject_id = None

        fhir_resources = process(results_payload_dict, self.args)
        self.assertEquals(len(fhir_resources), 5)
        subject = fhir_resources[0]
        self.assertDictEqual(subject, {
            'resourceType': 'Patient',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': 'project1'
                    }
                ]
            },
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
