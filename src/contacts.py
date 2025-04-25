"""
Author: Rafael Lemos - rafaellemos42@gmail.com
Date: 12/08/2024

License: MIT License
"""

from math import dist
from numpy import dot, arccos, degrees
from numpy.linalg import norm

from src.classes import Contact
from src.distances import distances
import src.conditions as conditions


def contact_detection(protein, context):
    """
    Detects contacts between atoms in a given protein.

    Args:
        protein (Protein): The protein object containing chain, residue and atom objects.

    Returns:
        list: A list of Contact objects representing the detected contacts.
    """
    region, interface, custom_distances, epsilon = context.region, context.interface, context.custom_distances, context.epsilon
    residues = list(protein.get_residues())
    contacts = []
    interface_res = set()
    max_ca_distance = 20.47 # 0.01 higher than the Arg-Arg pair
    
    count_contacts = {
        "hydrogen_bond":["HB",0],
        "hydrophobic":["HY",0],
        "attractive":["AT",0],
        "repulsive":["RE",0],
        "salt_bridge":["SB",0],
        "disulfide_bond":["DS",0],
        "stacking":["AS",0],
        "polar-apolar":["PA",0],
        "pos-apolar":["PosA",0],
        "neg-apolar":["NegA",0]
    }
    
    total_strength = 0
    
    contact_strength = {
        'hydrophobic': 0.6,
        'stacking': 1.5,
        'hydrogen_bond': 2.6,
        'attractive': 10.0,
        'repulsive': 10.0,
        'salt_bridge': 10.0,
        'disulfide_bond': 85.0,
        'polar-apolar': 0,
        'pos-apolar': 0,
        'neg-apolar': 0
    }
    
    residue_types = {
        'polar': ['C', 'H', 'N', 'Q', 'S', 'T', 'Y', 'W'],
        'apolar': ['A', 'F', 'G', 'I', 'V', 'M', 'P', 'L'],
        'charged': ['E', 'D', 'K', 'R']
    }
    
    residue_types_ic = {
        'polar': ['N', 'Q', 'S', 'T'],
        'apolar': ['A', 'C', 'G', 'F', 'I', 'M', 'L', 'P', 'W', 'V', 'Y'],
        'charged': ['E', 'D', 'H', 'K', 'R']
    }
    
    prodigy_count = {
        'charged-charged':0, 
        'charged-polar':0, 
        'charged-apolar':0,
        'polar-polar':0, 
        'polar-apolar':0,
        'apolar-apolar':0
    }
    prodigy = set()

    categories = custom_distances if custom_distances else conditions.categories
    if epsilon > 0:
        max_ca_distance += epsilon
        updated_distances = {key: value + epsilon for key, value in distances.items()}
    else:
        updated_distances = distances
        
    # if interface != 0:
    #     chains = interface.split(":")
    #     print(chains)
        
    for i, residue1 in enumerate(residues[1:]):
        for _, residue2 in enumerate(residues[i+1:], start=i+1):
            
            if residue1.resnum == residue2.resnum and residue1.chain.id == residue2.chain.id: # ignores same residue
                continue
            
            if region and (residue1.resnum not in region or residue2.resnum not in region):
                continue
            # if len(residue1.atoms) > 1 and len(residue2.atoms) > 1:
            #     ca1, ca2 = residue1.atoms[1], residue2.atoms[1] # alpha carbons

            #     distance_ca = dist((ca1.x, ca1.y, ca1.z), (ca2.x, ca2.y, ca2.z))
                
            #     # filter distant residues (static value then specific values)
            #     if distance_ca > max_ca_distance:
            #         continue
            #     else:
            #         key = ''.join(sorted((residue1.resname, residue2.resname)))
            #         if distance_ca > (updated_distances[key] + epsilon):
            #             continue
            # else:
            #     continue              
            
            # CHECKING FOR AROMATIC STACKINGS
            if residue1.ring and residue2.ring:
                ring1, ring2 = residue1.atoms[-1], residue2.atoms[-1] # RNG atoms
                if interface and ring1.entity == ring2.entity:
                    continue
                
                distance = dist((ring1.x, ring1.y, ring1.z), (ring2.x, ring2.y, ring2.z))
                angle = calc_angle(residue1.normal_vector, residue2.normal_vector)
                
                aromatic_range = categories['aromatic']
                if aromatic_range[0] <= distance <= aromatic_range[1]:
                    if (160 <= angle < 180) or (0 <= angle < 20):
                        stack_type = "-parallel"
                    elif (80 <= angle < 100):
                        stack_type = "-perpendicular"
                    else:
                        stack_type = "-other"

                    contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, ring1.atomname, 
                                    protein.id, residue2.chain.id, residue2.resnum, residue2.resname, ring2.atomname, 
                                    float(f"{distance:.2f}"), "stacking"+stack_type, ring1, ring2)
                    
                    count_contacts['stacking'][1] += 1
                    
                    contacts.append(contact)
                    
            for atom1 in residue1.atoms:
                for atom2 in residue2.atoms:

                    if interface and atom1.entity == atom2.entity:
                        continue
                    
                    name1 = f"{atom1.residue.resname}:{atom1.atomname}" # matches the pattern from conditions dictionary
                    name2 = f"{atom2.residue.resname}:{atom2.atomname}"

                    if name1 in conditions.contact_types and name2 in conditions.contact_types: # excludes the RNG atom and any different other
                        
                        distance = dist((atom1.x, atom1.y, atom1.z), (atom2.x, atom2.y, atom2.z))
                        if distance <= (6 + epsilon): # max distance for contacts

                            for contact_type, distance_range in categories.items():

                                if contact_type == 'hydrogen_bond' and (abs(residue2.resnum - residue1.resnum) <= 3): # skips alpha-helix for h-bonds
                                    continue
                                elif contact_type in ['polar-apolar', 'pos-apolar', 'neg-apolar'] and not interface:
                                    continue
                                
                                if distance_range[0] <= distance <= distance_range[1]: # fits the range
                                    if conditions.contact_conditions[contact_type](name1, name2): # fits the type of contact
                                        
                                        contact = Contact(protein.id, residue1.chain.id, residue1.resnum, residue1.resname, atom1.atomname, 
                                                        protein.id, residue2.chain.id, residue2.resnum, residue2.resname, atom2.atomname, 
                                                        float(f"{distance:.2f}"), contact_type, atom1, atom2)
                                        contacts.append(contact)
                                        interface_res.add(f"{residue1.chain.id},{residue1.resnum},{residue1.resname}")
                                        
                                        count_contacts[contact_type][1] += 1
                                        
                                        if interface: 
                                            total_strength += contact_strength[contact_type]
                        
                        # PRODIGY BLOCK
                        res1 = f"{residue1.chain.id}:{residue1.resnum}{residue1.resname}"
                        res2 = f"{residue2.chain.id}:{residue2.resnum}{residue2.resname}"
                        
                        if distance <= 5.5 and (res1, res2) not in prodigy:
                            # residue_group = next((k for k, v in residue_types_ic.items() if residue1.resname in v), None)
                            # print(residue1.resname, residue_group)
                            # print(res1, res2)
                            
                            if (residue1.resname in residue_types_ic['charged'] and residue2.resname in residue_types_ic['apolar']) or \
                            (residue1.resname in residue_types_ic['apolar'] and residue2.resname in residue_types_ic['charged']):
                                prodigy_count['charged-apolar'] += 1

                            elif (residue1.resname in residue_types_ic['polar'] and residue2.resname in residue_types_ic['apolar']) or \
                                (residue1.resname in residue_types_ic['apolar'] and residue2.resname in residue_types_ic['polar']):
                                prodigy_count['polar-apolar'] += 1

                            elif residue1.resname in residue_types_ic['charged'] and residue2.resname in residue_types_ic['charged']:
                                prodigy_count['charged-charged'] += 1

                            elif (residue1.resname in residue_types_ic['charged'] and residue2.resname in residue_types_ic['polar']) or \
                                (residue1.resname in residue_types_ic['polar'] and residue2.resname in residue_types_ic['charged']):
                                prodigy_count['charged-polar'] += 1

                            elif residue1.resname in residue_types_ic['polar'] and residue2.resname in residue_types_ic['polar']:
                                prodigy_count['polar-polar'] += 1

                            elif residue1.resname in residue_types_ic['apolar'] and residue2.resname in residue_types_ic['apolar']:
                                prodigy_count['apolar-apolar'] += 1
                            prodigy.add((res1, res2))

    values = prodigy_count.values()
    #print("id, charged-apolar, polar-apolar, charged-charged, charged-polar, polar-polar, apolar-apolar")
    print(f"{protein.id},"+",".join(str(v) for v in values))
                    
    if total_strength == 0:
        total_strength = None

    return contacts, interface_res, total_strength, count_contacts


def show_contacts(contacts):
    """
    Formats and summarizes contact information to be outputted to a file. Only works with the -o flag.

    Args:
        contacts (list): A list of Contact objects of a given protein.

    Returns:
        str: A formatted string summarizing the contact information.
    """
    
    output = []
    
    output.append("Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type")
    for contact in contacts:
        output.append(contact.print_text())
        
    return "\n".join(output) # returns as a string to be written directly into the file


def calc_angle(vector1, vector2):
    """
    Calculates the angle between two ring vectors of aromatic residues

    Args:
        vector1 (tuple): The first vector (x, y, z).
        vector2 (tuple): The second vector (x, y, z).

    Returns:
        float: The angle between the vectors in degrees.
    """
    
    dot_product = dot(vector1, vector2)
    magnitude_product = norm(vector1) * norm(vector2) # normalizes the dot product
    angle = arccos(dot_product / magnitude_product) # angle in radians   
    
    return degrees(angle)
