def initplot():
    fig = plt.figure(figsize=(plotparam["bbox"][0]/plotparam["dpi"], plotparam["bbox"][1]/plotparam["dpi"]), dpi=plotparam["dpi"])
    plt.axes().set_aspect('equal')
    plt.axes().set_xmargin(0.01)
    plt.axes().set_ymargin(0.01)
    plt.axes().set_axis_off()
    global zorder
    zorder = 0
    return fig


def nxdraw(G, networktype, map_center = False, nnids = False, objs = ["nodes", "edges"]):
    """Take an igraph graph G and draw it.
    The order of objects (nodes/edges) in objs determines zorder.
    """
    global zorder
    for obj in objs:
        if obj == "edges":
            drawfunc = "nx.draw_networkx_edges"
            plotparam = plotparam_edges
        elif obj == "nodes":
            drawfunc = "nx.draw_networkx_nodes"
            plotparam = plotparam_nodes
        G_nx = G.to_networkx()
            
        if nnids is not False: # Restrict to nnids node ids
            nnids_nx = [k for k,v in dict(G_nx.nodes(data=True)).items() if v['id'] in nnids]
            G_nx = G_nx.subgraph(nnids_nx)

        pos_transformed, map_center = project_nxpos(G_nx, map_center)
        coll = eval(drawfunc)(G_nx, pos_transformed, **plotparam[networktype])
        zorder += 1
        coll.set_zorder(zorder) # Draw edges on top of nodes. See: https://networkx.org/documentation/stable/_modules/networkx/drawing/nx_pylab.html#draw_networkx
    return map_center


def common_entries(*dcts):
    """Like zip() but for dicts.
    See: https://stackoverflow.com/questions/16458340/python-equivalent-of-zip-for-dictionaries
    """
    if not dcts:
        return
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def project_nxpos(G, map_center = False):
    """Take a spatial nx network G and projects its GPS coordinates to local azimuthal.
    Returns transformed positions, as used by nx.draw()
    """
    lats = nx.get_node_attributes(G, 'x')
    lons = nx.get_node_attributes(G, 'y')
    pos = {nid:(lat,-lon) for (nid,lat,lon) in common_entries(lats,lons)}
    if map_center:
        loncenter = map_center[0]
        latcenter = map_center[1]
    else:
        loncenter = np.mean(list(lats.values()))
        latcenter = -1* np.mean(list(lons.values()))
    local_azimuthal_projection = "+proj=aeqd +R=6371000 +units=m +lat_0={} +lon_0={}".format(latcenter, loncenter)
    # Use transformer: https://gis.stackexchange.com/questions/127427/transforming-shapely-polygon-and-multipolygon-objects
    wgs84_to_aeqd = pyproj.Transformer.from_proj(
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),
        pyproj.Proj(local_azimuthal_projection))
    pos_transformed = {nid:list(ops.transform(wgs84_to_aeqd.transform, Point(latlon)).coords)[0] for nid, latlon in pos.items()}
    return pos_transformed, (loncenter,latcenter)


def project_pos(lats, lons, map_center = False):
    """Project GPS coordinates to local azimuthal.
    """
    pos = [(lat,-lon) for lat,lon in zip(lats,lons)]
    if map_center:
        loncenter = map_center[0]
        latcenter = map_center[1]
    else:
        loncenter = np.mean(list(lats.values()))
        latcenter = -1* np.mean(list(lons.values()))
    local_azimuthal_projection = "+proj=aeqd +R=6371000 +units=m +lat_0={} +lon_0={}".format(latcenter, loncenter)
    # Use transformer: https://gis.stackexchange.com/questions/127427/transforming-shapely-polygon-and-multipolygon-objects
    wgs84_to_aeqd = pyproj.Transformer.from_proj(
        pyproj.Proj("+proj=longlat +datum=WGS84 +no_defs"),
        pyproj.Proj(local_azimuthal_projection))
    pos_transformed = [(ops.transform(wgs84_to_aeqd.transform, Point(latlon)).coords)[0] for latlon in pos]
    return pos_transformed, (loncenter,latcenter)


def round_coordinates(G, r = 7):
    for v in G.vs:
        G.vs[v.index]["x"] = round(G.vs[v.index]["x"], r)
        G.vs[v.index]["y"] = round(G.vs[v.index]["y"], r)

def mirror_y(G):
    for v in G.vs:
        y = G.vs[v.index]["y"]
        G.vs[v.index]["y"] = -y
    

def osm_to_ig(node, edge):
    """ Turns a node and edge dataframe into an igraph Graph.
    """
    
    G = ig.Graph(directed = False)

    x_coords = node['x'].tolist() 
    y_coords = node['y'].tolist()
    ids = node['osmid'].tolist()
    coords = []

    for i in range(len(x_coords)):
        G.add_vertex(x = x_coords[i], y = y_coords[i], id = ids[i])
        coords.append((x_coords[i], y_coords[i]))

    id_dict = dict(zip(G.vs['id'], np.arange(0, G.vcount()).tolist()))
    coords_dict = dict(zip(np.arange(0, G.vcount()).tolist(), coords))

    edge_list = []
    for i in range(len(edge)):
        edge_list.append([id_dict.get(edge['u'][i]), id_dict.get(edge['v'][i])])
        
    G.add_edges(edge_list)
    G.simplify()
    new_edges = G.get_edgelist()
    
    distances_list = []
    for i in range(len(new_edges)):
        distances_list.append(haversine(coords_dict.get(new_edges[i][0]), coords_dict.get(new_edges[i][1])))

    G.es()['weight'] = distances_list
    
    return G


def delete_overlaps(G_res, G_orig, verbose = False):
    """Deletes inplace all overlaps of G_res with G_orig (from G_res)
    based on node ids. In other words: G_res -= G_orig
    """
    del_edges = []
    for e in list(G_res.es):
        try:
            n1_id = e.source_vertex["id"]
            n2_id = e.target_vertex["id"]
            # If there is already an edge in the original network, delete it
            n1_index = G_orig.vs.find(id = n1_id).index
            n2_index = G_orig.vs.find(id = n2_id).index
            if G_orig.are_connected(n1_index, n2_index):
                del_edges.append(e.index)
        except:
            pass
    G_res.delete_edges(del_edges)
    # Remove isolated nodes
    isolated_nodes = G_res.vs.select(_degree_eq=0)
    G_res.delete_vertices(isolated_nodes)
    if verbose: print("Removed " + str(len(del_edges)) + " overlapping edges and " + str(len(isolated_nodes)) + " nodes.")

def intersect_igraphs(G1, G2):
    """Generates the graph intersection of igraph graphs G1 and G2, copying also link and node attributes.
    """
    # Ginter = G1.__and__(G2) # This does not work with attributes.
    if G1.ecount() > G2.ecount(): # Iterate through edges of the smaller graph
        G1, G2 = G2, G1
    inter_nodes = set()
    inter_edges = []
    inter_edge_attributes = {}
    inter_node_attributes = {}
    edge_attribute_name_list = G2.edge_attributes()
    node_attribute_name_list = G2.vertex_attributes()
    for edge_attribute_name in edge_attribute_name_list:
        inter_edge_attributes[edge_attribute_name] = []
    for node_attribute_name in node_attribute_name_list:
        inter_node_attributes[node_attribute_name] = []
    for e in list(G1.es):
        n1_id = e.source_vertex["id"]
        n2_id = e.target_vertex["id"]
        try:
            n1_index = G2.vs.find(id = n1_id).index
            n2_index = G2.vs.find(id = n2_id).index
        except ValueError:
            continue
        if G2.are_connected(n1_index, n2_index):
            inter_edges.append((n1_index, n2_index))
            inter_nodes.add(n1_index)
            inter_nodes.add(n2_index)
            edge_attributes = e.attributes()
            for edge_attribute_name in edge_attribute_name_list:
                inter_edge_attributes[edge_attribute_name].append(edge_attributes[edge_attribute_name])

    # map nodeids to first len(inter_nodes) integers
    idmap = {n_index:i for n_index,i in zip(inter_nodes, range(len(inter_nodes)))}

    G_inter = ig.Graph()
    G_inter.add_vertices(len(inter_nodes))
    G_inter.add_edges([(idmap[e[0]], idmap[e[1]]) for e in inter_edges])
    for edge_attribute_name in edge_attribute_name_list:
        G_inter.es[edge_attribute_name] = inter_edge_attributes[edge_attribute_name]

    for n_index in idmap.keys():
        v = G2.vs[n_index]
        node_attributes = v.attributes()
        for node_attribute_name in node_attribute_name_list:
            inter_node_attributes[node_attribute_name].append(node_attributes[node_attribute_name])
    for node_attribute_name in node_attribute_name_list:
        G_inter.vs[node_attribute_name] = inter_node_attributes[node_attribute_name]

    return G_inter


def ox_to_csv(G, p, placeid, parameterid, postfix = "", compress = True, verbose = True):
    if "crs" not in G.graph:
        G.graph["crs"] = 'epsg:4326' # needed for OSMNX's graph_to_gdfs in utils_graph.py
    try:
        node, edge = ox.graph_to_gdfs(G)
    except ValueError:
        node, edge = gpd.GeoDataFrame(), gpd.GeoDataFrame()
    prefix = placeid + '_' + parameterid + postfix

    node.to_csv(p + prefix + '_nodes.csv', index = False)
    if compress: compress_file(p, prefix + '_nodes')
 
    edge.to_csv(p + prefix + '_edges.csv', index = False)
    if compress: compress_file(p, prefix + '_edges')

    if verbose: print(placeid + ": Successfully wrote graph " + parameterid + postfix)


def check_extract_zip(p, prefix):
    """ Check if a zip file prefix+'_nodes.zip' and + prefix+'_edges.zip'
    is available at path p. If so extract it and return True, otherwise False.
    If you call this function, remember to clean up (i.e. delete the unzipped files)
    after you are done like this:

    if compress:
        os.remove(p + prefix + '_nodes.csv')
        os.remove(p + prefix + '_edges.csv')
    """

    try: # Use zip files if available
        with zipfile.ZipFile(p + prefix + '_nodes.zip', 'r') as zfile:
            zfile.extract(prefix + '_nodes.csv', p)
        with zipfile.ZipFile(p + prefix + '_edges.zip', 'r') as zfile:
            zfile.extract(prefix + '_edges.csv', p)
        return True
    except:
        return False

def compress_file(p, f, filetype = ".csv", delete_uncompressed = True):
    with zipfile.ZipFile(p + f + ".zip", 'w', zipfile.ZIP_DEFLATED) as zfile:
        zfile.write(p + f + filetype, f + filetype)
    if delete_uncompressed: os.remove(p + f + filetype)


def csv_to_ox(p, placeid, parameterid):
    """ Load a networkx graph from _edges.csv and _nodes.csv
    The edge file must have attributes u,v,osmid
    The node file must have attributes y,x,osmid
    Only these attributes are loaded, and edge lengths are calculated.
    """
    prefix = placeid + '_' + parameterid
    compress = check_extract_zip(p, prefix)
    
    with open(p + prefix + '_edges.csv', 'r') as f:
        header = f.readline().strip().split(",")

        lines = []
        for line in csv.reader(f, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
            line_list = [c for c in line]
            osmid = str(eval(line_list[header.index("osmid")])[0]) if isinstance(eval(line_list[header.index("osmid")]), list) else line_list[header.index("osmid")] # If this is a list due to multiedges, just load the first osmid
            line_string = "" + line_list[header.index("u")] + " "+ line_list[header.index("v")] + " " + osmid
            lines.append(line_string)
        G = nx.parse_edgelist(lines, nodetype = int, data = (("osmid", int),), create_using = nx.MultiDiGraph) # MultiDiGraph is necessary for OSMNX, for example for get_undirected(G) in utils_graph.py
    with open(p + prefix + '_nodes.csv', 'r') as f:
        header = f.readline().strip().split(",")
        values_x = {}
        values_y = {}
        for line in csv.reader(f, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
            line_list = [c for c in line]
            osmid = int(line_list[header.index("osmid")])
            values_x[osmid] = float(line_list[header.index("x")])
            values_y[osmid] = float(line_list[header.index("y")])

        nx.set_node_attributes(G, values_x, "x")
        nx.set_node_attributes(G, values_y, "y")

    edge_lengths = {}
    for e in G.edges(keys=True):
        edge_lengths[e] = haversine((values_x[e[0]], values_y[e[0]]), (values_x[e[1]], values_y[e[1]]))
    nx.set_edge_attributes(G, edge_lengths, "length")

    if compress:
        os.remove(p + prefix + '_nodes.csv')
        os.remove(p + prefix + '_edges.csv')
    return G


def csv_to_ig(p, placeid, parameterid):
    prefix = placeid + '_' + parameterid
    compress = check_extract_zip(p, prefix)
    empty = False

    try:
        n = pd.read_csv(p + prefix + '_nodes.csv')
        e = pd.read_csv(p + prefix + '_edges.csv')
    except:
        empty = True
    if compress:
        os.remove(p + prefix + '_nodes.csv')
        os.remove(p + prefix + '_edges.csv')
    if empty:
        return ig.Graph(directed = False)
    G = osm_to_ig(n, e)
    round_coordinates(G)
    mirror_y(G)
    return G


print("Loaded functions.\n")
