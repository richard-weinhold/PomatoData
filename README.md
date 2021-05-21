# pomato_data

This repository contains code for the creation of data for the POMATO 2030 Update.

# Data Sources: 
	Nodes, Lines: GridKit [done]
		Zubau: Manuell [Todo, aber pyPSA hat das eigentlich]
	Demand: ENTSO-E (zonal) und GDP-Disaggregation [done]
	Conventional Plants (OPSD, 2018) [todo: Manuelles decomminssioning]
	Renewables:
		Capacity: 
			Onshore/PV: Anymod Zonal, regionalisiert durch NUTS3 Potential [todo]
			Offshore: Anmod, Node-mapping Manuell
			Other: Teil von Bestand, Zubau Anymod Regionalisierung ?  
		Availability: Atlite NUTS3 [done], Offshore [FFE]

	Speicher: Teil von Conv, Zubau Manuell

	Sektorenkopplung (Teil von Anymod):
		Wäreme: ? 
		PtX: ? 

		Überlegung: Ist das flexibilität? 











