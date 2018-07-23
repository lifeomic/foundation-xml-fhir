#!/usr/bin/env python
import argparse
import json
import uuid
import xmltodict


def read_xml(xml_file):
    with open(xml_file) as fd:
        return xmltodict.parse(fd.read())


def save_json(fhir_resources, out_file):
    with open(out_file, 'wb') as fd:
        json.dump(fhir_resources, fd, indent=4)


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
    return specimen, specimen_id


def process(results_payload_dict, project_id, subject_id, out_file):
    fhir_resources = []
    specimen = create_specimen(
        results_payload_dict, project_id, subject_id)
    fhir_resources.append(specimen)
    save_json(fhir_resources, out_file)


def main():
    parser = argparse.ArgumentParser(
        prog='foundation-xml-fhir', description='Converts FoundationOne XML reports into FHIR resources.')
    parser.add_argument('-x, --xml', dest='xml_file',
                        required=True, help='Path to the XML file')
    parser.add_argument('-p, --project', dest='project_id', required=True,
                        help='The ID of the project to link the resources to')
    parser.add_argument('-s, --subject', dest='subject_id', required=True,
                        help='The ID of the subject/patient to link the resources to')
    parser.add_argument('-o, --output', dest='out_file',
                        required=True, help='Path to write the FHIR JSON resources')

    args = parser.parse_args()

    xml_dict = read_xml(args.xml_file)
    process(xml_dict['rr:ResultsReport']['rr:ResultsPayload'],
            args.project_id, args.subject_id, args.out_file)


if __name__ == "__main__":
    main()
