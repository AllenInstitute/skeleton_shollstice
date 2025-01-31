import os
import numpy as np
import pandas as pd

SWC_COLUMNS = (
    "id",
    "type",
    "x",
    "y",
    "z",
    "radius",
    "parent",
)
COLUMN_CASTS = {"id": int, "parent": int, "type": int}



class Shollster:

    def __init__(self, neuron, shell_thickness = 5):
        """
        
        Class for performing sholl analysis


        Args:
            neuron (neuron_morphology.morphology or path to an swc file)
            shell_thickness (float) size of shell thickness for sholl analysis
            
        """
        self.neuron = neuron
        self.shell_thickness = shell_thickness
        self.entire_morph_df = load_neuron(self.neuron)
        process_morph_df(self.entire_morph_df, self.shell_thickness)
        self.analysis_morph_df = self.entire_morph_df
    
    def run_sholl(self, sholl_method , normalize = False, compatment_types = [2, 3, 4,] ):
        """
        conduct sholl analysis
        
        Args:

            sholl_method (str): ['intersections', 'nodes', 'length', 'branches' ]
            
            shell_thickness (float): the size of concentric spheres 
            
            normalize (bool): wether to divide values by the area/volume of the sholl ring
            
        Returns:
            quantified sholl analysis  
        """
        # make sure soma is in there
        compatment_types = compatment_types + [1]

        self.analysis_morph_df = self.entire_morph_df[self.entire_morph_df['type'].isin(compatment_types)]
        sorted_shells = [i for  i in sorted( set(self.entire_morph_df.shell.values), key=lambda x: int(x.split("_")[0]))]

        if sholl_method == "intersections":
            
            return _intersection_sholl_analysis(self.analysis_morph_df, sorted_shells, self.shell_thickness,  normalize)
        
        elif sholl_method == 'length':
            return _total_length_sholl_analysis(self.analysis_morph_df, sorted_shells,  normalize)

        elif sholl_method == "nodes":
            return _num_nodes_sholl_analysis(self.analysis_morph_df, sorted_shells,  normalize)
            
        elif sholl_method == "branches":
            return _branches_sholl_analysis(self.analysis_morph_df, sorted_shells,  normalize)

        else:
            supported_methods = ['intersections', 'nodes', 'length', 'branches' ]
            raise ValueError("Unsupported sholl method provided: {}. Supported methods are:\n{}".format(sholl_method,supported_methods))

def load_neuron(neuron):
    
    if type(neuron)==str:
        if not os.path.exists(neuron):
            raise ValueError(f"Path to skeleton does not exist:\n{neuron}")
        if not neuron.endswith(".swc"):
            raise ValueError(f"Can only read .swc files at this time")
                
        morph_df = pd.read_csv(neuron, names=SWC_COLUMNS, comment="#", sep=" ", index_col=False)
        for key,typ in COLUMN_CASTS.items():
            morph_df[key] = morph_df[key].astype(typ)
            
        return morph_df
    
    elif hasattr(neuron, 'nodes'):
        morph_df = pd.DataFrame(neuron.nodes())
        return morph_df
    else:
        neuron_type = type(neuron)
        raise ValueError(f"expecting either neuron_morphology.Morphology object or path to an swc file for `neuron`. Received:\n{neuron_type}")


def process_morph_df(morph_df, shell_thickness):
    
    "identify node types and distances from soma"
    
    child_counts = morph_df['parent'].value_counts().to_dict()
    morph_df["num_children"] = morph_df.id.map(child_counts).fillna(0)

    # Identify all node IDs
    all_nodes = set(morph_df['id'])
    parent_nodes = set(child_counts.keys())
    tip_nodes = all_nodes - parent_nodes

    # Assign node types
    morph_df['node_type'] = morph_df['id'].map(lambda node: 
        "branch" if child_counts.get(node, 0) > 1 else 
        "tip" if node in tip_nodes else 
        "reducible"
    )
    
    node_positions = morph_df.set_index('id')[['x', 'y', 'z']].to_dict('index')

    # Compute distance from the soma (node with parent -1)
    soma_node_id = morph_df.loc[ (morph_df['parent'] == -1) & (morph_df['type'] == 1), 'id'].values[0]
    soma_coords = node_positions[soma_node_id]

    morph_df['distance_from_soma'] = morph_df['id'].map(lambda node: 
        np.linalg.norm(np.array(list(node_positions[node].values())) - np.array(list(soma_coords.values())))
    )
    
    # Compute segment lengths (Euclidean distance to parent)
    def segment_length(node, parent):
        if parent == -1:
            return 0  # Root node has no segment length
        return np.linalg.norm(np.array(list(node_positions[node].values())) - np.array(list(node_positions[parent].values())))

    morph_df['segment_length'] = morph_df.apply(lambda row: segment_length(row['id'], row['parent']), axis=1)

    # # Bin segment lengths into shells
    morph_df['shell'] = (morph_df['distance_from_soma'] // shell_thickness).astype(int)
    morph_df['shell']  = morph_df['shell'].apply(lambda x: f"{x}_{x * shell_thickness}_{(x + 1) * shell_thickness}")


def _branches_sholl_analysis(morph_df, sorted_shells, normalize):
    
    morph_df = morph_df[morph_df['node_type']=='branch']
    num_nodes_per_shell = pd.DataFrame(morph_df.groupby('shell').size(), columns=['num_branches'])
    num_nodes_per_shell = num_nodes_per_shell.reindex(sorted_shells, fill_value=0)
    # num_nodes_per_shell = num_nodes_per_shell.loc[sorted_shells]
    
    return num_nodes_per_shell


def _num_nodes_sholl_analysis(morph_df, sorted_shells, normalize):
    
    num_nodes_per_shell = pd.DataFrame(morph_df.groupby('shell').size(), columns=['num_nodes'])
    num_nodes_per_shell = num_nodes_per_shell.reindex(sorted_shells, fill_value=0)
    # num_nodes_per_shell = num_nodes_per_shell.loc[sorted_shells]
    return num_nodes_per_shell


def _total_length_sholl_analysis(morph_df, sorted_shells, normalize):
    
    shell_lengths = pd.DataFrame(morph_df.groupby('shell')['segment_length'].sum())
    shell_lengths = shell_lengths.reindex(sorted_shells, fill_value=0)
    # shell_lengths = shell_lengths.loc[sorted_shells]
    
    return shell_lengths    

def _intersection_sholl_analysis(morph_df, sorted_shells, shell_thickness, normalize):
    """
    Count number of intersections in each ring

    Args:
        morph_df (pd.DataFrame): processed dataframe 
        sorted_shells (list): list of shell names (format: 'index_shell0_shell0+1') in ascending order from soma
        shell_thickness (float): thickness of the shells used
        normalize (bool): normalize shell by area/volume

    Returns:
        pd.DataFrame: _description_
    """
    intersection_count = 0
    intersection_records = {}

    for shell in sorted_shells:
        shell_df = morph_df[morph_df['shell']==shell]
        branch_df = shell_df[shell_df['node_type']=='branch']
        tip_df = shell_df[shell_df['node_type']=='tip']
        
        if not branch_df.empty:
            
            if shell == f'0_0_{shell_thickness}':
                intersection_count += branch_df.num_children.sum() 
            else:
                intersection_count += (branch_df.num_children.sum() -  branch_df.shape[0])
            
        intersection_count -= tip_df.shape[0]
        
        intersection_records[shell] = intersection_count

    res = pd.DataFrame(intersection_records.items(),columns=['shell','num_intersections'])
    res=res.set_index('shell')
    return res


