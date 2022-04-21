Hallo,
your job on rarefan.evolbio.mpg.de with ID {{job.run_id}} is complete.

You can access the results at http://rarefan.evolbio.mpg.de/results?run_id={{job.run_id}}.

## Job Summary
{% if job.overall_status == "complete with warnings" %}
There have been warnings or errors during the postprocessing of your
results. Please inspect the output data and logfile (out/rarefan.log)
carefully.
{% endif %}


### RAYTs
{% if job.stages['rarefan']['results']['data_sanity']['rayts'] == 1 %}. **There was a problem with the RAYT analysis output data. Please check your results carefully.** 

{% endif %} We discovered {{job.stages['rarefan']['results']['counts']['rayts']}} RAYTs across all submitted genomes using tblastn with {{job.setup["query_rayt"]}} at an e-value threshold of {{job.setup["e_value_cutoff"]}}.


### Seed sequences
{% if job.stages['rarefan']['results']['data_sanity'] == 1 %} **There was a problem with the seed sequence output data. Please check your results carefully.**

{% endif %} There are {{job.stages['rarefan']['results']['counts']['nmers']}} {{job.setup['nmer_length']}}bp long sequences in the reference genome that occur more frequently than {{job.setup['min_nmer_occurrence']}} times.


### REPINs
{% if job.stages['rarefan']['results']['data_sanity']['repins'] == 1 %}. **There was a problem with the REP(IN) analysis output data. Please check your results carefully.**

{% endif %} {% if job.setup['analyse_repins'] in ['1', 1, 'y', True] %} We detected {{job.stages['rarefan']['results']['counts']['repins']}} REPINs in the reference genome. {% else %} We detected {{job.stages['rarefan']['results']['counts']['repins']}} REPs in the reference genome. {% endif %}

In case of problems, please reply to this email.

Thank you for using RAREFAN.

http://rarefan.evolbio.mpg.de
