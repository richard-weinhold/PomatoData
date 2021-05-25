# pomato_data

This repository contains code for the creation of data for the POMATO 2030 Update.

# Data Sources: 
	Nodes, Lines: GridKit [done]
		Zubau: 
			- DC Lines (TYNDP), source GridKit https://github.com/PyPSA/pypsa-eur/blob/master/data/links_tyndp.csv 
				- As file add_dclines.csv which adds new lines to the network. 

	Demand: OPSD (ENTSO-E) (zonal) und GDP-Disaggregation [done]

	Conventional Plants (OPSD, 2018) [todo: Manuelles decomminssioning]

	Renewables:
		Capacity: 
			Onshore/PV: Anymod Zonal, regionalisiert durch NUTS3 Potential [todo]
			Offshore: Anmod, Node-mapping Manuell
			Other: Teil von Bestand, Zubau Anymod Regionalisierung ?  
		Availability: Onshore/PV Atlite NUTS3 [done], Offshore Atlite EEZ. 


	Speicher: Teil von Conv, Zubau Manuell
		- Inflow Timeseries aus Atlite hydro. 
			- TODO: Umrechnung von m^3 pro stunde zu MWh

	Sektorenkopplung (Teil von Anymod):
		Wäreme: ? 
		PtX: ? 

		Überlegung: Ist das flexibilität? 











