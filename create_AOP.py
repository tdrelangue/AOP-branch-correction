from SPARQLWrapper import SPARQLWrapper, JSON
import json
from visualize_AOP import print


def aop_dump(aop_id):    # endpoint sparql    
    """
    Queries the AOP-Wiki SPARQL endpoint for information related to a specific AOP.

    Parameters:
    - aop_id (int or str): Must be an integer or a string containing only digits.
      Represents the AOP identifier (e.g., 131).

    Returns:
    - dict: JSON-style dictionary following the SPARQLWrapper format.
      Contains bindings for Molecular Initiating Events (MIE), Key Events (KEs),
      Key Event Relationships (KERs), Adverse Outcomes (AO), and their associated labels,
      titles, genes, and identifiers.

    Notes:
    - If the query fails, the exception is printed.
    - The endpoint used is: https://aopwiki.rdf.bigcat-bioinformatics.org/sparql
    """
    sparql = SPARQLWrapper("https://aopwiki.rdf.bigcat-bioinformatics.org/sparql")    
    if isinstance(aop_id, int) or (isinstance(aop_id, str) and aop_id.isdigit()):        
        sparql.setReturnFormat(JSON)        
        sparql.setQuery(f"""SELECT DISTINCT ?AOP ?MIE ?KE_up ?KE_dwn ?AO ?ker_genes ?aop_id ?aop_label ?ke_id ?ke_label ?ke_title ?ke_genes ?ke_dwn_label ?ke_dwn_id ?ke_dwn_title ?ke_dwn_genes 
                        WHERE{{  BIND(aop:{aop_id} AS ?aop_id)  ?AOP a aopo:AdverseOutcomePathway ;       
                        dc:identifier ?aop_id;       
                        rdfs:label ?aop_label .  OPTIONAL {{ ?AOP aopo:has_molecular_initiating_event ?MIE . }}  OPTIONAL {{ ?AOP aopo:has_adverse_outcome ?AO . }}  ?AOP aopo:has_key_event ?KE_up .  ?KE_up dc:identifier ?ke_id ;
                        rdfs:label ?ke_label ;         
                        dc:title ?ke_title .  OPTIONAL {{ ?KE_up edam:data_1025 ?ke_genes . }}  OPTIONAL {{    ?AOP aopo:has_key_event_relationship ?KER .    
                        ?KER a aopo:KeyEventRelationship ;         
                        aopo:has_upstream_key_event ?KE_up ;         
                        aopo:has_downstream_key_event ?KE_dwn .    #OPTIONAL {{ ?KER edam:data_1025 ?ker_genes . }}    
                        ?KE_dwn dc:identifier ?ke_dwn_id ;            rdfs:label ?ke_dwn_label ;            dc:title ?ke_dwn_title .    OPTIONAL {{ ?KE_dwn edam:data_1025 ?ke_dwn_genes . }}  }}}}""")        
        try:            
            ret = sparql.query()            
            json_format = ret.convert()            
            return json_format        
        except Exception as e:            
            print(e)


def add_AOP_variable_by_keid(AOP,new_var_dict, new_var_name="AC50"):
    """
    Adds a new variable to each node in the AOP graph based on matching KE_id values.

    Parameters:
    - AOP (dict): Dictionary representing the AOP graph. Each key is a node name (e.g., 'KE0'),
      and each value is a dict with keys like 'KE_id', 'name', 'connections', and 'genes'.
    - new_var_dict (dict): Dictionary mapping KE_id (str) -> value (any type, e.g., float or str).
    - new_var_name (str): Name of the variable to add (e.g., "AC50", "probability").

    Returns:
    - dict: The updated AOP dictionary with the new variable added to corresponding nodes.
    """
    for node in AOP:
        for key in new_var_dict:
            if AOP[node]["KE_id"]==key:
                AOP[node][new_var_name]= new_var_dict[key]
    return AOP


def build_AOP(result, MIES, AOS, ordered_nodes):
    """
    Constructs a structured AOP graph dictionary from SPARQL query results.

    Parameters:
    - result (dict): SPARQL result in JSON format containing KE metadata.
    - MIES (list of str): List of KE_id strings that are Molecular Initiating Events.
    - AOS (list of str): List of KE_id strings that are Adverse Outcomes.
    - ordered_nodes (list of str): Ordered list of KE_id strings based on causal flow.

    Returns:
    - dict: A dictionary representing the AOP graph. Keys are formatted node names (e.g., 'MIE0', 'KE1', 'AO0'),
      and values are dicts containing:
        - 'name' (str): KE title,
        - 'connections' (list): Empty list, to be filled later,
        - 'genes' (list of str): Gene identifiers associated with the KE (if available),
        - 'KE_id' (str): Original KE identifier.
    """

    AOP = {}
    for i, node in enumerate(ordered_nodes):
        key=f"KE{i}"
        if node in AOS:
            key = f"AO{AOS.index(node)}"
        if node in MIES:
            key = f"MIE{MIES.index(node)}"
        ke_title = ""
        ke_id = node
        ke_genes = []
        genesseen = set()
        for binding in result["results"]["bindings"]:
            label = binding["ke_label"]["value"]
            title = binding["ke_title"]["value"]
            try:
                gene = binding["ke_genes"]["value"].split("/")[-1]
                if label == f"KE {node}":
                    ke_title = title
                    if gene not in genesseen:
                        ke_genes.append(gene)
                        genesseen.add(gene)
            except:
                #who cares it is not used yet anyway
                continue
        AOP[key] = {
            "name": ke_title,
            "connections": [], 
            "genes": ke_genes,
            "KE_id": ke_id
        }
        
    return AOP


def get_nodes_in_apparition_order(result, connections):
    """
    Determines the KE_id node order based on traversal starting from the MIE node.

    Parameters:
    - result (dict): SPARQL JSON result; used only to extract the first MIE KE_id.
    - connections (list of tuples): Each tuple contains two KE_id strings (upstream, downstream).

    Returns:
    - list of str: Ordered list of KE_id values representing the causal flow of events in the AOP.
    """

    temp_connections = [connection for connection in connections]
    ordered_nodes = [result["results"]["bindings"][0]["MIE"]["value"].split("/")[-1]]
    i=0
    seen = set()  # This will keep track of items we've already seen
    while temp_connections:
        previous_node = ordered_nodes[i]
        connections_left=[]
        for connection in temp_connections:
            if connection[0] == previous_node:
                if connection[1] not in seen:
                    ordered_nodes.append(connection[1])
                    seen.add(connection[1])
            else :
                connections_left.append(connection)

        # Let's not forget the useless ones that don't have MIE as an ancestor
        temp_connections = connections_left
        connections_left=[]
        for connection in temp_connections:
            if connection[1] == previous_node:
                if connection[0] not in seen:
                    ordered_nodes.append(connection[0])
                    seen.add(connection[0])
            else :
                connections_left.append(connection)
        temp_connections = connections_left
        i+=1
    return ordered_nodes


def find_AO_and_MIE(result):
    """
    Extracts unique lists of Molecular Initiating Events (MIEs) and Adverse Outcomes (AOs) from the SPARQL result.

    Parameters:
    - result (dict): SPARQL result as a JSON-style dictionary.

    Returns:
    - tuple:
        - MIES (list of str): KE_ids that are marked as MIEs.
        - AOS (list of str): KE_ids that are marked as AOs.
    """

    seen_AO = set()
    seen_MIE = set()
    MIES =[]
    AOS = []
    for binding in result["results"]["bindings"]:
        node = binding["MIE"]["value"].split("/")[-1]
        if node not in seen_MIE:
            MIES.append(node)
            seen_MIE.add(node)
        
        node = binding["AO"]["value"].split("/")[-1]
        if node not in seen_AO:
            AOS.append(node)
            seen_AO.add(node)
    return MIES,AOS


def collect_connections_in_AOP(result):
    """
    Extracts all upstream-downstream relationships (connections) between KEs from the SPARQL result.

    Parameters:
    - result (dict): SPARQL query result in JSON format.

    Returns:
    - list of tuple(str, str): Each tuple represents a connection (KE_up_id, KE_down_id).
      Ensures all pairs are unique. Tuples are always of length 2.
    """

    connections = []
    seen = set()  # This will keep track of items we've already seen
    for binding in result["results"]["bindings"]:
        if "KE_dwn" in binding:
            conn = (binding['KE_up']["value"].split("/")[-1],
                    binding['KE_dwn']["value"].split("/")[-1]
                    )
            
            if conn not in seen:
                connections.append(conn)
                seen.add(conn)
    return connections


def add_connections_to_AOP(connections, AOP):
    """
    Adds downstream connections to each node in the AOP structure.

    Parameters:
    - connections (list of tuple(str, str)): Each tuple is (KE_up_id, KE_down_id).
    - AOP (dict): Dictionary representing the AOP, with nodes keyed by name and containing at least a 'KE_id'.

    Operation:
    - For each connection, the function searches for both KE_id matches in the AOP,
      and appends the destination node name to the 'connections' list of the source node.

    Returns:
    - None (the AOP dict is modified in-place).
    """

    for connection in connections:
        for node in AOP:
            if connection[1]==AOP[node]["KE_id"]:
                connection = (connection[0], node)
        for node in AOP:
            if connection[0]==AOP[node]["KE_id"]:
                AOP[node]["connections"].append(connection[1])


def create_AOP_from_scratch(aop_id, manualKEEdges=None):
    """
    Main wrapper function to create a fully structured AOP graph from an AOP ID.

    Parameters:
    - aop_id (int or str): A valid AOP identifier.

    Steps:
    - Queries the SPARQL endpoint for the given AOP.
    - Collects KER connections.
    - Identifies MIEs and AOs.
    - Determines node order from MIE.
    - Builds node metadata (name, genes, KE_id).
    - Fills in connections between nodes.

    Returns:
    - dict: Final AOP graph with all nodes and their connections.
    """

    # Step 1: Retrieve raw AOP data
    result = aop_dump(aop_id)

    # Step 2: Extract connection 
    connections = collect_connections_in_AOP(result)
    if manualKEEdges:
        print(f"Using manual KE edges: {manualKEEdges}")
        for connection in manualKEEdges:
            connections.append(connection)

    # Step 3: collect and order the nodes in order of apparition
    MIES, AOS = find_AO_and_MIE(result)
    ordered_nodes = get_nodes_in_apparition_order(result, connections)

    # Step 4: Build base AOP structure
    AOP = build_AOP(result, MIES, AOS, ordered_nodes)

    # Step 5: Add connections (manual if provided, else default)
    add_connections_to_AOP(connections, AOP)
    return AOP

def add_proba_by_keid(AOP, proba):
    """
    Adds posterior probability values to each node in the AOP based on KE_id.

    Parameters:
    - AOP (dict): The AOP graph structure.
    - proba (dict): Dictionary mapping KE_id -> probability value (e.g., float between 0 and 1).

    Returns:
    - None (the AOP dict is modified in-place).
    """

    add_AOP_variable_by_keid(AOP=AOP, new_var_dict=proba, new_var_name="P(prior|event)")

if __name__ == "__main__":
    aop_id = "372"  # Replace with a valid AOP ID
    AOP = create_AOP_from_scratch(aop_id=aop_id,manualKEEdges=[('26', '1614'), ('1614', '286'), ('286', '1616'), ('1616', '1839')])


    print(AOP)