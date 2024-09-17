"""
    This is the main script to convert a mechanism to SSH-aerosol v2.0 format.
"""

import os

from parameters import (
    reaction_type,
    reaction_file,
    species_file,
    species_type,
    chem_id,
    output_dir,
    tag_fake,
    soap_file,
)
from data_stream import read_chem_sets, to_ssh_sets


if __name__ == "__main__":

    print("Start conversion...")
    # Read the mechanism
    if reaction_type == "GECKO":
        print("Reading GECKO mechanism from folder:", reaction_file)
    elif reaction_type == "FACSMILE":
        print("Reading MCM mechanism from files:", reaction_file, species_file)
    elif reaction_type == "SSH":
        print("Reading SSH mechanism from files:", reaction_file, species_file)
    else:
        raise ValueError("Unknown reaction type:", reaction_type)
    reactions, species = read_chem_sets(
        reaction_file, species_file, reaction_type, species_type, aero_vfile=soap_file
    )

    # Modification of the mechanism is done here if needed
    # ...

    # Output the mechanism in SSH-aerosol v2.0 format
    print(
        f"Outputting mechanism {chem_id} to folder:", os.path.join(output_dir, chem_id)
    )
    to_ssh_sets(output_dir, chem_id, reactions, species, tag_fake=tag_fake)
    print("Conversion is done!")
