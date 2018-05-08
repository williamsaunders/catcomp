William Saunders, May 7, 2018

This outlines the steps to create Planet Nine detection efficiency statistics for 
catalog comparison. 

This process creates many P9s to simulate, simulates them, determines observations on DES
determines which observations are realized at which magnitude, the number of requisite 
detections, and whether those detections constitute a discovery. 

Steps: 

1. sphere_points.py 
   - Create a set of starting sphere points in lat, lon (RA, dec)
   - Divide those into 10 separate, random files 
   - Outputs .csv files into sphere_points/

2. p9population.py
   - Simulate one chunk of those sphere points
   - Can qsub into X nodes using qsub_p9.py
   - Simulates the orbits and determines "observations" of P9 on individual CCDs
   - Outputs .fits tables into P9simulation_results

3. recombine_p9.py
   - After qsubbing p9population.py, recombines separate chunks into one file

4. combine_tile_folio.py
   - Combines tiles into one .fits table for each zone
   - Required to combine because exposure area is smaller than tile area
   - Can qsub with qsub-combine_tiles.py
   - Fast (~minutes) but uses a LOT of memory

5. detection_efficiency_folio.py 
   - Calculates the point-source detection efficiency for a zone
   - Can feed it multiple zones to single-thread
   - Can qsub with qsub=det_eff.py
   - Outputs .fits tables into zone_efficiencies
   - Use the efficiency outputs to determine whether an observation is "realized" 
     at a specific P9 magnitude

6. recombine_det-prob.py
   - Takes each zone efficiency and combines into one table unique for each exposure
   - Overlapping exposures are chosen on the basis of # of coadd objects (keeps exposure 
     with more)

7. discovery.py
   - Total() function calculates total number of observed P9s
   - Discovery() function perfors discovery. Takes a P9 magnitude (99=bright limit)
     and a number of detections that constitutes a discovery
   - Discovery() looks up point-source efficiency for the exposure to determine if realized

8. P9discovery.py 
   - Uses functions in discovery.py
   - Writes all discovery results into pickle dump files labeled .fits, but they are not .fits

9. P9discovery_plot.py 
   - Uses results from P9discovery.py 
   - Ccreates discovery plot result and histogram of # of detections

10. m50hist.py 
   - Creates histogram of m50 values for each exposure

