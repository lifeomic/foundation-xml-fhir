from unittest import TestCase
from src.convert import process

results_payload_dict = {
    'variant-report': {
        'samples': {
            'sample': {
                '@name': 'sample1'
            }
        }
    },
    'FinalReport': {
        'PMI': {
            'LastName': 'doe',
            'FirstName': 'jane',
            'MRN': '1234',
            'Gender': 'female',
            'DOB': '1970-01-01',
            'CollDate': '2000-01-01',
            'SubmittedDiagnosis': 'Cancer'
        },
        'Sample': {
            'TestType': 'Test 1'
        }
    }
}

class ConvertTest(TestCase):
    def test_convert_with_subject(self):
        self.maxDiff = 1200
        fhir_resources = process(results_payload_dict, 'project1', 'subject1')
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
            'variant': []
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
            'specimen': {
                'display': 'sample1',
                'reference': 'Specimen/{}'.format(specimen['id'])
            },
            'result': [],
            'id': report['id']
        })

    def test_convert_with_no_subject(self):
        fhir_resources = process(results_payload_dict, 'project1', None)
        self.assertEquals(len(fhir_resources), 4)
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
        self.assertEquals(specimen['subject']['reference'], 'Patient/{}'.format(subject['id']))
        sequence = fhir_resources[2]
        self.assertEquals(sequence['patient']['reference'], 'Patient/{}'.format(subject['id']))
        report = fhir_resources[3]
        self.assertEquals(report['subject']['reference'], 'Patient/{}'.format(subject['id']))

