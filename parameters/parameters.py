plotparam = {
			"bbox": (1280,1280),
			"dpi": 96
			}

plotparam_nodes = {
			"carall": {"node_color": '#918d87', "node_size": 6*4},
			"bike_onstreet": {"node_color": '#2222ff', "node_size": 3*1.5},
			"bike_offstreet": {"node_color": '#00991f', "node_size": 3*1.5},
			"bike_sharedstreet": {"node_color": '#ff0000', "node_size": 6*4},
			"car30": {"node_color": '#bbbbbb', "node_size": 6*4}
			}

plotparam_edges = {
			"carall": {"width": 6, "edge_color": '#918d87'},
			"bike_onstreet": {"width": 4.5, "edge_color": '#2222ff'},
			"bike_offstreet": {"width": 4.5, "edge_color": '#00991f'},
			"bike_sharedstreet": {"width": 6, "edge_color": '#ff0000'},
			"car30": {"width": 6, "edge_color": '#bbbbbb'}
			}

osmnxparameters = {'car30': {'network_type':'drive', 'custom_filter':'["maxspeed"~"^30$|^20$|^15$|^10$|^5$|^20 mph|^15 mph|^10 mph|^5 mph"]', 'retain_all': True},
                   'carall': {'network_type':'drive', 'custom_filter': None, 'retain_all': True},
                   'bike_cyclewaytrack': {'network_type':'bike', 'custom_filter':'["cycleway"~"track"]', 'retain_all': True},
                   'bike_highwaycycleway': {'network_type':'bike', 'custom_filter':'["highway"~"cycleway"]', 'retain_all': True},
                   'bike_designatedpath': {'network_type':'all', 'custom_filter':'["highway"~"path"]["bicycle"~"designated"]', 'retain_all': True},
                   'bike_cyclewayrighttrack': {'network_type':'bike', 'custom_filter':'["cycleway:right"~"track"]', 'retain_all': True},
                   'bike_cyclewaylefttrack': {'network_type':'bike', 'custom_filter':'["cycleway:left"~"track"]', 'retain_all': True},
                   'bike_cyclestreet': {'network_type':'bike', 'custom_filter':'["cyclestreet"]', 'retain_all': True},
                   'bike_bicycleroad': {'network_type':'bike', 'custom_filter':'["bicycle_road"]', 'retain_all': True},
                   'bike_livingstreet': {'network_type':'bike', 'custom_filter':'["highway"~"living_street"]', 'retain_all': True}
                  }  
# https://wiki.openstreetmap.org/wiki/Key:cycleway#Cycle_tracks
# https://wiki.openstreetmap.org/wiki/Tag:highway=path#Usage_as_a_universal_tag
# https://wiki.openstreetmap.org/wiki/Tag:highway%3Dliving_street
# https://wiki.openstreetmap.org/wiki/Key:cyclestreet


print("Loaded parameters.\n")
