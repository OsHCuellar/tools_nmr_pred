import io
import json
import logging

import flask
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError
from compute.gaussianCalc import gaussPoints, importSDF#, pandasDF

blueprint = flask.Blueprint("compute", __name__, url_prefix="/compute")

logger = logging.getLogger("tools-app")


@blueprint.route("/process_structure/", methods=["POST"])
def process_structure():
    """Example view to process a crystal structure."""

    # check if the post request has the file part,
    # otherwise redirect to first page
    if "structurefile" not in flask.request.files:
        # This will redirect the user to the selection page,
        # that is called `input_data` in tools-barebone
        return flask.redirect(flask.url_for("input_data"))

    # Get structure, file format, file content, and form data
    # (needed for additional information, e.g. cell in the case
    # of a XYZ file)
    structurefile = flask.request.files["structurefile"]
    fileformat = flask.request.form.get("fileformat", "unknown")
    filecontent = structurefile.read().decode("utf-8")
    fileobject = io.StringIO(str(filecontent))
    form_data = dict(flask.request.form)

    # Use
    try:
        structure_tuple = get_structure_tuple(
            fileobject, fileformat, extra_data=form_data
        )
    except UnknownFormatError:
        # You can use the flask.flash functionality to send a message
        # back to the structure selection page; this
        # will be shown in a red box on the top
        flask.flash("Unknown format '{}'".format(fileformat))
        return flask.redirect(flask.url_for("input_data"))
    except Exception:
        # Let's deal properly with any exception, to avoid to get a 500 error.
        # Feel free to do better error management here,
        # or to pass additional information via flask.flash
        flask.flash(
            "I tried my best, but I wasn't able to load your "
            "file in format '{}'...".format(fileformat)
        )
        return flask.redirect(flask.url_for("input_data"))
    # If we are here, the file was retrieved.
    # It will contain a tuple of length three, with:
    # - the 3x3 unit cell (in angstrom)
    # - a Nx3 list of atomic coordinates (in angstrom)
    # - a list of integer atomic numbers of length N

    # As an example, we just create a string representation of the JSON
    # and send it back to the user, to be rendered in a form
    data_for_template = {
        "structure_json": json.dumps(
            {
                "cell": structure_tuple[0],
                "atoms": structure_tuple[1],
                "numbers": structure_tuple[2],
            },
            indent=2,
            sort_keys=True,
        )
    }
    return flask.render_template(
        "user_templates/tools_example.html", **data_for_template
    )

@blueprint.route("/prediction/", methods=["GET", "POST"])
def prediction():
    """Processing the endpoint /compute/process_example_value/.

    This is called by the additional HTML form added in user_templates/additional_selection.html.

    Note that this is just for testing purposes: in reality you want to
    use `flask.render_template` as in the example above, rather than
    just returning a string, as you will need to pass the
    full HTML headers."""
    if flask.request.method == "POST":

        flask.flash("The prediction backend is not implemented yet.")
        return flask.redirect(flask.url_for("input_data"))

    return flask.redirect(flask.url_for("input_data"))

@blueprint.route("/chemical_diversity/", methods=["GET", "POST"])
def process_structure_example():
    """Processing the endpoint /compute/process_example_value/.

    This is called by the additional HTML form added in user_templates/additional_selection.html.

    Note that this is just for testing purposes: in reality you want to
    use `flask.render_template` as in the example above, rather than
    just returning a string, as you will need to pass the
    full HTML headers."""
    if flask.request.method == "POST":

        minX=-5
        maxX=5
        npoints=50
        a=1 #height
        b=0 #mean value
        c=1 #standard deviation

        #gaussPoints(minX, maxX, npoints, a, b, c)
        graphX, graphY, graphZ, molFormula, molSmiles, molWeight = importSDF()
        #df = pandasDF()

        #graphX, graphY = gaussPoints(minX,maxX,npoints,a,b,c)
        #graphZ = graphX + graphY
        data_for_template = {
                "graphX": str(graphX),
                "graphY": str(graphY),
                "minX": str(minX),
                "maxX": str(maxX),
                "XplusY": str(graphZ),
                "molFormula": str(molFormula),
                "molSmiles": str(molSmiles),
                "molWeight": str(molWeight),
                }
        #return str(graphX) + " ## " + str(graphY) + " ## " + str(graphZ)

        # This is the `value` of the option tag, not the text shown to the user
        value = flask.request.form.get("examplevalue", "UNKNOWN")
        return flask.render_template("user_templates/result_example.html", **data_for_template)
        return "This was a POST request, with value <pre>{}</pre>".format(value)
    return "This was a GET request"
