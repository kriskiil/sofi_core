from __future__ import absolute_import

import connexion


from . import encoder


def script_runs_within_container():
    with open("/proc/1/cgroup", "r") as cgroup_file:
        return "docker" in cgroup_file.read()

def create_app():
    options = {"swagger_ui": True}
    app = connexion.App(__name__, specification_dir='./openapi/', options=options)
    app.app.json_encoder = encoder.JSONEncoder
    app.add_api('openapi.yaml',
                arguments={'title': 'SAP'},
                pythonic_params=True)
    return app.app


if __name__ == '__main__':
    create_app().run(port=8080)
