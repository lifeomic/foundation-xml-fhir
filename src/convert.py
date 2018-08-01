#!/usr/bin/env python
import argparse
import json
import logging
import uuid
import xmltodict
import os
import gzip
import shutil
from utils import parse_hgvs


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


def create_observation(fasta, genes, project_id, subject_id, specimen_id, specimen_name, sequence_id):
    def create(variant_dict):
        observation_id = str(uuid.uuid4())
        position_value = variant_dict['@position']
        region, position = position_value.split(':')
        transcript = variant_dict['@transcript']
        cds_effect = variant_dict['@cds-effect'].replace('&gt;', '>')
        variant_name = '{}:c.{}'.format(transcript, cds_effect)
        chrom, offset, ref, alt = parse_hgvs(variant_name, fasta, genes)
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
                    }
                ]
            },
            'status': variant_dict['@status'],
            'subject': {
                'reference': 'Patient/{}'.format(subject_id)
            },
            'specimen': {
                'display': specimen_name,
                'reference': 'Specimen/{}'.format(specimen_id)
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
                    'url': 'http://hl7.org/fhir/StructureDefinition/observation-geneticsSequence',
                    'valueReference': {
                        'reference': 'Sequence/{}'.format(sequence_id)
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
        return observation
    return create


def create_report(results_payload_dict, project_id, subject_id, specimen_id, specimen_name, file_url=None):
    report_id = str(uuid.uuid4())

    report = {
        'resourceType': 'DiagnosticReport',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
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
            'test': results_payload_dict['FinalReport']['Sample']['TestType']
        },
        'issued': results_payload_dict['FinalReport']['PMI']['CollDate'],
        'subject': {
            'reference': 'Patient/{}'.format(subject_id)
        },
        'specimen': [{
            'display': specimen_name,
            'reference': 'Specimen/{}'.format(specimen_id)
        }],
        'result': [],
        'id': report_id
    }

    if file_url is not None:
        report['presentedForm'] = {
            'url': file_url,
            'contentType': 'application/pdf',
            'title': results_payload_dict['FinalReport']['Sample']['TestType']
        }

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


def create_specimen(results_payload_dict, project_id, subject_id):
    specimen_name = results_payload_dict['variant-report']['samples']['sample']['@name']
    specimen_id = str(uuid.uuid4())

    specimen = {
        'resourceType': 'Specimen',
        'meta': {
            'tag': [
                {
                    'system': 'http://lifeomic.com/fhir/dataset',
                    'code': project_id
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


def process(results_payload_dict, args):
    fhir_resources = []
    subject_id = args.subject_id

    if subject_id is None:
        subject, subject_id = create_subject(
            results_payload_dict, args.project_id)
        fhir_resources.append(subject)

    specimen, specimen_id, specimen_name = create_specimen(
        results_payload_dict, args.project_id, subject_id)
    sequence, sequence_id = create_sequence(
        args.project_id, subject_id, specimen_id, specimen_name)
    report = create_report(results_payload_dict, args.project_id,
                           subject_id, specimen_id, specimen_name, args.file_url)

    observations = []
    if ('short-variants' in results_payload_dict['variant-report'].keys() and
            'short-variant' in results_payload_dict['variant-report']['short-variants'].keys()):
        observations = list(map(create_observation(args.fasta, args.genes, args.project_id, subject_id, specimen_id, specimen_name, sequence_id),
                                results_payload_dict['variant-report']['short-variants']['short-variant']))

    report['result'] = [
        {'reference': 'Observation/{}'.format(x['id'])} for x in observations]
    sequence['variant'] = [
        {'reference': 'Observation/{}'.format(x['id'])} for x in observations]

    fhir_resources.append(specimen)
    fhir_resources.append(sequence)
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


if __name__ == '__main__':
    main()
