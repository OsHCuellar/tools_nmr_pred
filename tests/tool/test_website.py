import os

import pytest

from selenium.webdriver.support.ui import Select
from urllib.parse import urlparse

# Note: you need to build docker and start it with
# ./admin-tools/build-docker.sh && ./admin-tools/run-docker.sh
# Then, you can run the tests (e.g. with )
TEST_URL = "http://localhost:8091"
STRUCTURE_EXAMPLES_PATH = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "xsf_structure_examples"
)


@pytest.mark.nondestructive
def test_binput_data_page(selenium):
    """Check the page that is shown by default.

    Check that it contains the structure selector from tools-barebone."""
    selenium.get(TEST_URL)

    assert "tools-example" in selenium.title.lower()
    format_selector = selenium.find_element_by_id("fileformatSelect")

    # This is not a complete list, but at least these should be present
    expected_importer_names = set(
        [
            "Quantum ESPRESSO input [parser: qe-tools]",
            "VASP POSCAR [parser: ase]",
            "XCrySDen (.xsf) [parser: ase]",
            "CIF File (.cif) [parser: pymatgen]",
            "XYZ File (.xyz) [parser: ase]",
        ]
    )

    # If the difference is not empty, at least one of the expected importer names is not there!
    assert not expected_importer_names.difference(
        option.text for option in format_selector.find_elements_by_tag_name("option")
    )

    # Check the presence of a string in the source code
    assert "YES-val" in selenium.page_source


def get_file_examples():
    """Get all valid files from the STRUCTURE_EXAMPLES_PATH and returns
    a list of filename relative paths."""
    retval = []
    for filename in os.listdir(STRUCTURE_EXAMPLES_PATH):
        if filename.endswith("~") or filename.startswith("."):
            continue
        retval.append(filename)

    return retval


def submit_xsf_structure(selenium, file_abspath):
    """Given a selenium driver, submit a file."""
    # Load file
    file_upload = selenium.find_element_by_name("structurefile")
    file_upload.send_keys(file_abspath)

    # Select format
    format_selector = selenium.find_element_by_id("fileformatSelect")
    Select(format_selector).select_by_value("xsf-ase")

    # Submit form
    # selenium.find_element_by_xpath("//input[@value='Calculate my structure']").click()
    selenium.find_element_by_xpath(
        "//form[@action='compute/process_structure/']"
    ).submit()


@pytest.mark.nondestructive
@pytest.mark.parametrize("file_relpath", get_file_examples())
def test_send_structure(selenium, file_relpath):
    """Test submitting various files."""
    selenium.get(TEST_URL)

    # Load file
    file_abspath = os.path.join(STRUCTURE_EXAMPLES_PATH, file_relpath)
    submit_xsf_structure(selenium, file_abspath)

    # We should have been redirected back to /
    assert urlparse(selenium.current_url).path == "/compute/process_structure/"

    assert "Successfully parsed structure tuple" in selenium.page_source


@pytest.mark.nondestructive
@pytest.mark.parametrize(
    "example_value, example_value_readable_string",
    [["YES-val", "Yes"], ["NO-val", "No"], ["MAYBE-val", "Maybe"],],
)
def test_send_example_value(selenium, example_value, example_value_readable_string):
    """Test submitting various test values."""
    selenium.get(TEST_URL)

    # Select format
    format_selector = selenium.find_element_by_name("examplevalue")
    Select(format_selector).select_by_visible_text(example_value_readable_string)

    # Submit form
    # selenium.find_element_by_xpath("//input[@value='Calculate my structure']").click()
    selenium.find_element_by_xpath(
        "//form[@action='compute/process_example_value/']"
    ).submit()

    # We should be on the right page (and no redirection should have occurred)
    assert urlparse(selenium.current_url).path == "/compute/process_example_value/"

    assert example_value in selenium.page_source
