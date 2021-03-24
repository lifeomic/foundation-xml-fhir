#!/usr/bin/env python
import argparse
import base64
import json
import logging
import uuid
import xmltodict
import os
import gzip
import datetime
import shutil
from utils import parse_hgvs, parse_splice
from subprocess import call


logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s')
logger = logging.getLogger(__name__)


def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def save_json(fhir_resources, out_file):
    with open(out_file, 'wb') as fd:
        json.dump(fhir_resources, fd, indent=4)


def unzip(zipped_file):
    unzipped_file = os.path.splitext(zipped_file)[0]
    logger.info('Unzipping %s to %s', zipped_file, unzipped_file)

    with gzip.open(zipped_file, "rb") as f_in, open(unzipped_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    logger.info('Unzipping completed')
    return unzipped_file


def create_microsatallite_observation(project_id, subject_id, specimen_id, effective_date, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())

        values = {
            'MSI-H': 'MSI-H',
            'MSI-L': 'MSI-L',
            'MSS': 'Stable',
            'unknown': 'unknown'
        }

        codes = {
            'MSI-H': 'LA26203-2',
            'MSI-L': 'LA26202-4',
            'MSS': 'LA14122-8',
            'unknown': 'unknown'
        }

        observation = {
            'resourceType': 'Observation',
            'effectiveDateTime': effective_date,
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': project_id
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'microsatellite-instability'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
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
            'status': 'final',
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
            'valueCodeableConcept': {
                'coding': [
                    {
                        'system': 'http://loinc.org',
                        'code': codes.get(variant_dict['@status'], 'unknown'),
                        'display': values.get(variant_dict['@status'], 'unknown')
                    }
                ]
            },
            'extension': [

            ],
            'id': observation_id
        }

        if specimen_id is not None:
            observation['specimen'] = {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
            }

        if sequence_id is not None:
             observation['extension'].append({
                'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                'valueReference': {
                    'reference': 'Sequence/{}'.format(sequence_id)
                }
            })
        return observation
    return create


def create_tumor_mutation_observation(project_id, subject_id, specimen_id, effective_date, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())

        codes = {
            'high': 'TMB-H',
            'intermediate': 'TMB-I',
            'low': 'TMB-L',
            'unknown': 'unknown'
        }

        observation = {
            'resourceType': 'Observation',
            'effectiveDateTime': effective_date,
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': project_id
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'tumor-mutation-burden'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
                    }
                ]
            },
            'code': {
                'coding': [
                {
                    'system': 'http://lifeomic.com/fhir/biomarker',
                    'code': 'TMB',
                    'display': 'Tumor Mutation Burden'
                }
                ]
            },
            'status': 'final',
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
             'extension': [

            ],
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
                                'code': codes.get(variant_dict['@status'], 'unknown'),
                                'display': variant_dict['@status']
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
                        'value': float(variant_dict['@score']),
                        'unit': variant_dict['@unit']
                    }
                }
            ],
            'id': observation_id
        }

        if specimen_id is not None:
            observation['specimen'] = {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
            }

        if sequence_id is not None:
            observation['extension'].append({
                'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                'valueReference': {
                    'reference': 'Sequence/{}'.format(sequence_id)
                }
            })
        return observation
    return create


def create_rearrangement_observation(project_id, subject_id, specimen_id, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())

        observation = {
            'resourceType': 'Observation',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': project_id
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'rearrangement'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
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
            },
            'status': 'final',
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': variant_dict['@status'],
                    'display': 'Foundation - {}'.format(variant_dict['@status'].title())
                    }
                ]
            },
            'extension': [
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://www.genenames.org',
                                'code': '1100',
                                'display': variant_dict.get('@targeted-gene', variant_dict.get('@target-gene'))
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48004-6',
                                'display': '{} {}'.format(variant_dict.get('@targeted-gene', variant_dict.get('@target-gene')), variant_dict['@type'].title())
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48002-0',
                                'display': 'somatic'
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48001-2',
                                'display': '{} {}'.format(variant_dict['@pos1'], variant_dict['@pos2'])
                            }
                        ]
                    }
                }
            ],
            'id': observation_id
        }

        if specimen_id is not None:
            observation['specimen'] = {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
            }

        if sequence_id is not None:
            observation['extension'].append({
                'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                'valueReference': {
                    'reference': 'Sequence/{}'.format(sequence_id)
                }
            })
        return observation
    return create


def create_copy_number_observation(project_id, subject_id, specimen_id, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())
        position_value = variant_dict['@position']
        region, position = position_value.split(':')

        observation = {
            'resourceType': 'Observation',
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': project_id
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                     {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'copy-number'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
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
            },
            'status': 'final',
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': variant_dict['@status'],
                    'display': 'Foundation - {}'.format(variant_dict['@status'].title())
                    }
                ]
            },
            'extension': [
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://www.genenames.org',
                                'code': '1100',
                                'display': variant_dict['@gene']
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48004-6',
                                'display': '{}: CN={}'.format(variant_dict['@type'].title(), variant_dict['@copy-number'])
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48002-0',
                                'display': 'somatic'
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48001-2',
                                'display': position
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAChromosome',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '47999-8',
                                'display': region
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsCopyNumberEvent',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://www.sequenceontology.org',
                                'code': 'SO:0001019',
                                'display': variant_dict['@type'].capitalize()
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeName',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48005-3',
                                'display': 'Exons {}'.format(variant_dict['@number-of-exons'])
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-copyNumber',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://lifeomic.com',
                                'code': 'copyNumber',
                                'display': variant_dict['@copy-number']
                            }
                        ]
                    }
                }
            ],
            'id': observation_id
        }

        if specimen_id is not None:
            observation['specimen'] = {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
            }

        if sequence_id is not None:
            observation['extension'].append({
                'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                'valueReference': {
                    'reference': 'Sequence/{}'.format(sequence_id)
                }
            })
        return observation
    return create


def hgvs_2_vcf (variant_name, genes, functional_effect, cds_effect, position_value, strand, fasta):
    if functional_effect in ['splice', 'frameshift', 'nonframeshift']:
      if cds_effect in ['2169_*27&gt;T','2169_*27>T']:
        cds_effect = '2169_*27del48'
      return parse_splice(cds_effect, position_value, strand, fasta)
    else:
        try:
            return parse_hgvs(variant_name, fasta, genes)
        except:
            return parse_splice(cds_effect, position_value, strand, fasta)


def create_observation(fasta, genes, project_id, subject_id, specimen_id, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())
        position_value = variant_dict['@position']
        functional_effect = variant_dict['@functional-effect']
        strand = variant_dict['@strand']
        region, position = position_value.split(':')
        transcript = variant_dict['@transcript']
        cds_effect = variant_dict['@cds-effect'].replace('&gt;', '>')
        variant_name = '{}:c.{}'.format(transcript, cds_effect)
        chrom, offset, ref, alt = hgvs_2_vcf(variant_name, genes, functional_effect, cds_effect, position_value, strand, fasta)
        variantReadCount = int(round(int(variant_dict['@depth']) * float(variant_dict['@allele-fraction'])))

        observation = {
            'resourceType': 'Observation',
            'identifier': [{
                'system': 'https://lifeomic.com/observation/genetic',
                'value': '{}:{}:{}:{}'.format(chrom, offset, ref, alt)
            }],
            'meta': {
                'tag': [
                    {
                        'system': 'http://lifeomic.com/fhir/dataset',
                        'code': project_id
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/source',
                        'code': 'LifeOmic Task Service'
                    },
                     {
                        'system': 'http://lifeomic.com/fhir/variant-type',
                        'code': 'short'
                    },
                    {
                        'system': 'http://lifeomic.com/fhir/report-source',
                        'code': 'Foundation'
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
            },
            'status': 'final',
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
            'valueCodeableConcept': {
                'coding': [
                    {
                    'system': 'http://foundationmedicine.com',
                    'code': variant_dict['@status'],
                    'display': 'Foundation - {}'.format(variant_dict['@status'].title())
                    }
                ]
            },
            'extension': [
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGene',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://www.genenames.org',
                                'code': '1100',
                                'display': variant_dict['@gene']
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsDNASequenceVariantName',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48004-6',
                                'display': variant_name
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeType',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://snomed.info/sct',
                                'code': 'LL380-7',
                                'display': variant_dict['@functional-effect']
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAminoAcidChangeName',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48005-3',
                                'display': 'p.{}'.format(variant_dict['@protein-effect'])
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsAllelicFrequency',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '81258-6',
                                'display': variant_dict['@allele-fraction']
                            }
                        ]
                    }
                },
                {
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsGenomicSourceClass',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48002-0',
                                'display': 'somatic'
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAPosition',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '48001-2',
                                'display': position
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsDNAChromosome',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '47999-8',
                                'display': region
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsTotalReadDepth',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '82121-5',
                                'display': variant_dict['@depth']
                            }
                        ]
                    }
                },
                {
                    "url": "http://lifeomic.com/fhir/StructureDefinition/observation-geneticsVariantReadCount",
                    "valueCodeableConcept": {
                        "coding": [
                            {
                                "system": "http://loinc.org",
                                "code": "82121-5",
                                "display": str(variantReadCount)
                            }
                        ]
                    }
                },
                {
                    'url': 'http://lifeomic.com/fhir/StructureDefinition/observation-geneticsTranscriptID',
                    'valueCodeableConcept': {
                        'coding': [
                            {
                                'system': 'http://loinc.org',
                                'code': '51958-7',
                                'display': variant_dict['@transcript']
                            }
                        ]
                    }
                }
            ],
            'id': observation_id
        }

        if specimen_id is not None:
            observation['specimen'] = {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
            }

        if sequence_id is not None:
            observation['extension'].append({
                'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                'valueReference': {
                    'reference': 'Sequence/{}'.format(sequence_id)
                }
            })
        return observation
    return create


def create_report(results_payload_dict, project_id, subject_id, specimen_id, specimen_name, effective_date,
                  file_url=None, sequence_id=None):
    report_id = str(uuid.uuid4())

    report = {
        'resourceType': 'DiagnosticReport',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
                },
                {
                    'system': 'http://lifeomic.com/fhir/source',
                    'code': 'LifeOmic Task Service'
                },
                {
                    'system': 'http://lifeomic.com/fhir/report-source',
                    'code': 'Foundation'
                }
            ]
        },
        'extension': [
            {
                'url': 'http://hl7.org/fhir/StructureDefinition/DiagnosticReport-geneticsAssessedCondition',
                'valueReference': {
                    'reference': results_payload_dict['FinalReport']['PMI']['SubmittedDiagnosis']
                }
            }
        ],
        'status': 'final',
        'code': {
            'text': results_payload_dict['FinalReport']['Sample']['TestType']
        },
        'effectiveDateTime': effective_date,
        'subject': {
            'reference': 'Patient/{}'.format(subject_id)
        },
        'result': [],
        'id': report_id
    }

    if specimen_id is not None:
        report['specimen'] = [{
            'display': specimen_name,
            'reference': 'Specimen/{}'.format(specimen_id)
        }]

    if file_url is not None:
        report['presentedForm'] = [{
            'url': file_url,
            'contentType': 'application/pdf',
            'title': results_payload_dict['FinalReport']['Sample']['TestType']
        }]

    if sequence_id is not None:
        report['extension'].append({
            'url': 'http://lifeomic.com/fhir/StructureDefinition/sequence-id',
            'valueString': sequence_id
        })

    return report


def create_subject(results_payload_dict, project_id):
    subject_id = str(uuid.uuid4())
    pmi_dict = results_payload_dict['FinalReport']['PMI']

    subject = {
        'resourceType': 'Patient',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
                },
                {
                    'system': 'http://lifeomic.com/fhir/source',
                    'code': 'LifeOmic Task Service'
                }
            ]
        },
        'name': [{
            'use': 'official',
            'family': pmi_dict['LastName'],
            'given': [pmi_dict['FirstName']]
        }],
        'identifier': [{
            'type': {
                'coding': [{
                    'system': 'http://hl7.org/fhir/v2/0203',
                    'code': 'MR'
                }]
            },
            'value': pmi_dict['MRN']
        }],
        'gender': pmi_dict['Gender'].lower(),
        'birthDate': pmi_dict['DOB'],
        'id': subject_id
    }
    return subject, subject_id


def create_sequence(project_id, subject_id, specimen_id, specimen_name):
    sequence_id = str(uuid.uuid4())

    sequence = {
        'resourceType': 'Sequence',
        'type': 'dna',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
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
        'patient': {
            'reference': 'Patient/{}'.format(subject_id)
        },
        'specimen': {
            'display': specimen_name,
            'reference': 'Specimen/{}'.format(specimen_id)
        },
        'referenceSeq': {
            'genomeBuild': 'GRCh37'
        },
        'id': sequence_id,
        'variant': []
    }
    return sequence, sequence_id


def get_specimen_name(results_payload_dict):
    specimen_name = None
    if isinstance(results_payload_dict['variant-report']['samples']['sample'], list):
        found = list(filter(lambda x:  x['@nucleic-acid-type'] == 'DNA', results_payload_dict['variant-report']['samples']['sample']))
        if len(found) > 0:
            specimen_name = found[0]['@name']
    else:
        specimen_name = results_payload_dict['variant-report']['samples']['sample']['@name']
    return specimen_name


def create_specimen(results_payload_dict, project_id, subject_id):
    specimen_name = get_specimen_name(results_payload_dict)
    specimen_id = str(uuid.uuid4())

    specimen = {
        'resourceType': 'Specimen',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
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
                'value': specimen_name
            }
        ],
        'subject': {
            'reference': 'Patient/{}'.format(subject_id)
        },
        'id': specimen_id
    }
    return specimen, specimen_id, specimen_name


def write_vcf(variants, specimen_name, fasta, genes, vcf_out_file):
    with open('./unsorted.vcf', 'w+') as vcf_file:
        vcf_file.write('##fileformat=VCFv4.2\n')
        vcf_file.write('##source=foundation-xml-fhir\n')
        vcf_file.write('##reference=file://{}\n'.format(fasta))
        vcf_file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
        vcf_file.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">\n')
        vcf_file.write('##INFO=<ID=VENDSIG,Number=1,Type=String,Description="Vendor Significance">\n')
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        vcf_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        vcf_file.write('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Number of reads harboring allele (in order specified by GT)">\n')
        vcf_file.write('##contig=<ID=chr1,length=248956422>\n')
        vcf_file.write('##contig=<ID=chr2,length=242193529>\n')
        vcf_file.write('##contig=<ID=chr3,length=198295559>\n')
        vcf_file.write('##contig=<ID=chr4,length=190214555>\n')
        vcf_file.write('##contig=<ID=chr5,length=181538259>\n')
        vcf_file.write('##contig=<ID=chr6,length=170805979>\n')
        vcf_file.write('##contig=<ID=chr7,length=159345973>\n')
        vcf_file.write('##contig=<ID=chr8,length=145138636>\n')
        vcf_file.write('##contig=<ID=chr9,length=138394717>\n')
        vcf_file.write('##contig=<ID=chr10,length=133797422>\n')
        vcf_file.write('##contig=<ID=chr11,length=135086622>\n')
        vcf_file.write('##contig=<ID=chr12,length=133275309>\n')
        vcf_file.write('##contig=<ID=chr13,length=114364328>\n')
        vcf_file.write('##contig=<ID=chr14,length=107043718>\n')
        vcf_file.write('##contig=<ID=chr15,length=101991189>\n')
        vcf_file.write('##contig=<ID=chr16,length=90338345>\n')
        vcf_file.write('##contig=<ID=chr17,length=83257441>\n')
        vcf_file.write('##contig=<ID=chr18,length=80373285>\n')
        vcf_file.write('##contig=<ID=chr19,length=58617616>\n')
        vcf_file.write('##contig=<ID=chr20,length=64444167>\n')
        vcf_file.write('##contig=<ID=chr21,length=46709983>\n')
        vcf_file.write('##contig=<ID=chr22,length=50818468>\n')
        vcf_file.write('##contig=<ID=chrX,length=156040895>\n')
        vcf_file.write('##contig=<ID=chrY,length=57227415>\n')
        vcf_file.write('##contig=<ID=chrM,length=16569>\n')
        vcf_file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(specimen_name))

        status = {
            'known': 'Pathogenic',
            'likley': 'Likely_pathogenic',
            'unknown': 'Uncertain_significance',
            'ambiguous': 'other'
        }

        for variant_dict in variants:
            vendsig = status.get(variant_dict.get('@status', 'unknown'))
            cds_effect = variant_dict['@cds-effect'].replace('&gt;', '>')
            transcript = variant_dict['@transcript']
            functional_effect = variant_dict['@functional-effect']
            strand = variant_dict['@strand']
            position_value = variant_dict['@position']
            dp = variant_dict['@depth']
            af = variant_dict['@allele-fraction']
            gt = '1/1' if float(variant_dict['@allele-fraction']) > 0.9 else '0/1'
            alt = int(round(int(dp) * float(af)))
            ref = int(dp) - alt
            ad = '{},{}'.format(ref, alt)
            variant_name = '{}:c.{}'.format(transcript, cds_effect)

            chrom, offset, ref, alt = hgvs_2_vcf(variant_name, genes, functional_effect, cds_effect, position_value, strand, fasta)
            vcf_file.write('{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={};AF={};VENDSIG={}\tGT:DP:AD\t{}:{}:{}\n'.format(chrom, offset, ref, alt, dp, af, vendsig, gt, dp, ad))


def process(results_payload_dict, args):
    fhir_resources = []
    subject_id = args.subject_id

    if subject_id is None:
        subject, subject_id = create_subject(
            results_payload_dict, args.project_id)
        fhir_resources.append(subject)

    specimen_name = None
    specimen_id = None
    sequence_id = None

    if (args.vcf_out_file is None):
        specimen, specimen_id, specimen_name = create_specimen(
            results_payload_dict, args.project_id, subject_id)
        sequence, sequence_id = create_sequence(
            args.project_id, subject_id, specimen_id, specimen_name)
        fhir_resources.append(specimen)
        fhir_resources.append(sequence)

    effective_date = datetime.datetime.now().isoformat()
    if 'CollDate' in results_payload_dict['FinalReport']['PMI'] and '#text' in results_payload_dict['FinalReport']['PMI']['CollDate']:
        effective_date = results_payload_dict['FinalReport']['PMI']['CollDate']['#text']
    elif 'CollDate' in results_payload_dict['FinalReport']['PMI']:
        effective_date = results_payload_dict['FinalReport']['PMI']['CollDate']

    report = create_report(results_payload_dict, args.project_id,
                           subject_id, specimen_id, specimen_name, effective_date, args.file_url, args.sequence_id)

    observations = []
    if ('short-variants' in results_payload_dict['variant-report'].keys()):

        variants = []

        if (results_payload_dict['variant-report']['short-variants'] is not None and
            'short-variant' in results_payload_dict['variant-report']['short-variants'].keys()):
            variants_dict = results_payload_dict['variant-report']['short-variants']['short-variant']
            variants = variants_dict if isinstance(variants_dict, list) else [variants_dict]

        if (args.vcf_out_file is not None):
            specimen_name = get_specimen_name(results_payload_dict)
            write_vcf(variants, specimen_name, args.fasta, args.genes, args.vcf_out_file)

        observations = list(map(create_observation(args.fasta, args.genes, args.project_id, subject_id, specimen_id, specimen_name, sequence_id or args.sequence_id),
                            variants))

    if ('copy-number-alterations' in results_payload_dict['variant-report'].keys()):

        cnvs = []

        if (results_payload_dict['variant-report']['copy-number-alterations'] is not None and
            'copy-number-alteration' in results_payload_dict['variant-report']['copy-number-alterations'].keys()):
            cnv_dict = results_payload_dict['variant-report']['copy-number-alterations']['copy-number-alteration']
            cnvs = cnv_dict if isinstance(cnv_dict, list) else [cnv_dict]

        observations.extend(list(map(create_copy_number_observation(args.project_id, subject_id, specimen_id, specimen_name, sequence_id or args.sequence_id),
                                cnvs)))

    if ('rearrangements' in results_payload_dict['variant-report'].keys()):

        rearrangements = []

        if (results_payload_dict['variant-report']['rearrangements'] is not None and
            'rearrangement' in results_payload_dict['variant-report']['rearrangements'].keys()):
            rearrangement_dict = results_payload_dict['variant-report']['rearrangements']['rearrangement']
            rearrangements = rearrangement_dict if isinstance(rearrangement_dict, list) else [rearrangement_dict]

        observations.extend(list(map(create_rearrangement_observation(args.project_id, subject_id, specimen_id, specimen_name, sequence_id or args.sequence_id),
                                rearrangements)))
    if ('biomarkers' in results_payload_dict['variant-report'].keys()):

        if (results_payload_dict['variant-report']['biomarkers'] is not None and
            'microsatellite-instability' in results_payload_dict['variant-report']['biomarkers'].keys()):
            microsatellite_dict = results_payload_dict['variant-report']['biomarkers']['microsatellite-instability']
            observations.append(create_microsatallite_observation(args.project_id, subject_id, specimen_id, effective_date, specimen_name, sequence_id or args.sequence_id)(microsatellite_dict))

        if (results_payload_dict['variant-report']['biomarkers'] is not None and
            'tumor-mutation-burden' in results_payload_dict['variant-report']['biomarkers'].keys()):
            tumor_dict = results_payload_dict['variant-report']['biomarkers']['tumor-mutation-burden']
            observations.append(create_tumor_mutation_observation(args.project_id, subject_id, specimen_id, effective_date, specimen_name, sequence_id or args.sequence_id)(tumor_dict))


    report['result'] = [
        {'reference': 'Observation/{}'.format(x['id'])} for x in observations]

    if (args.vcf_out_file is None):
        sequence['variant'] = [
            {'reference': 'Observation/{}'.format(x['id'])} for x in observations]

    fhir_resources.append(report)
    fhir_resources = fhir_resources + observations
    logger.info('Created %d FHIR resources', len(fhir_resources))

    return fhir_resources


def main():
    parser = argparse.ArgumentParser(
        prog='foundation-xml-fhir', description='Converts FoundationOne XML reports into FHIR resources.')
    parser.add_argument('-r, --reference', dest='fasta',
                        required=True, help='Path to reference genome')
    parser.add_argument('-g, --genes', dest='genes',
                        required=False, help='Path to genes file', default='/opt/app/refGene.hg19.txt')
    parser.add_argument('-x, --xml', dest='xml_file',
                        required=True, help='Path to the XML file')
    parser.add_argument('-p, --project', dest='project_id', required=True,
                        help='The ID of the project to link the resources to')
    parser.add_argument('-s, --subject', dest='subject_id', required=False,
                        help='The ID of the subject/patient to link the resources to')
    parser.add_argument('-o, --output', dest='out_file',
                        required=True, help='Path to write the FHIR JSON resources')
    parser.add_argument('-f, --file', dest='file_url',
                        required=False, help='The URL to the PDF Report in the PHC')
    parser.add_argument('-d, --pdf-output', dest='pdf_out_file',
                        required=False, help='Path to write the PDF file', default=None)
    parser.add_argument('-v, --vcf-output', dest='vcf_out_file',
                        required=False, help='Path to write the VCF file', default=None)
    parser.add_argument('-i, --sequence-id', dest='sequence_id',
                        required=False, help='The sequence id to add to the Diagnostic Report', default=None)

    args = parser.parse_args()
    logger.info('Converting XML to FHIR with args: %s',
                json.dumps(args.__dict__))

    # pyfaidx has a bug with bgzipped files.  Unzip the genome for now
    # https://github.com/mdshw5/pyfaidx/issues/125
    if (args.fasta.lower().endswith('.bgz') or
            args.fasta.lower().endswith('.gz')):
        args.fasta = unzip(args.fasta)

    xml_dict = read_xml(args.xml_file)
    fhir_resources = process(
        xml_dict['rr:ResultsReport']['rr:ResultsPayload'], args)
    save_json(fhir_resources, args.out_file)
    logger.info('Saved FHIR resources to %s', args.out_file)

    if args.pdf_out_file is not None and 'ReportPDF' in xml_dict['rr:ResultsReport']['rr:ResultsPayload']:
        pdf = base64.b64decode(xml_dict['rr:ResultsReport']['rr:ResultsPayload']['ReportPDF'])
        with open(args.pdf_out_file, "w") as pdf_file:
            pdf_file.write(pdf)
        logger.info('Saved PDF report to %s', args.pdf_out_file)

    if args.vcf_out_file is not None:
        call(['/opt/app/sort.sh', args.vcf_out_file])


if __name__ == '__main__':
    main()
