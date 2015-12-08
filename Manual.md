## `orthAgogue`, an advanced tool with an easy command line interface ##
`orthAgogue` is a tool for orthology estimation. The tool solves issues related to memory management, data reading speed, data storage and parallelization. An overview of these issues will appear later, a brief description of the tool can be found in a submitted paper. (For questions, contact the author at [oekseth@gmail.com](mailto:oekseth@gmail.com).)

Here we provide a list of Frequently Asked Questions. In order to make this FAQ both complete and relevant, your feedback is very much welcome. Therefore, if you have any questions, please forward them to [orthagogue-issue-tracker@googlegroups.com](mailto:orthagogue-issue-tracker@googlegroups.com).

**Q**: How did you manage to get a 100x speedup?
<br>
<b>A</b>: Yes, sound too good to be true, it was not a simple task. We focused on balancing effective data structures with parallelization to achieve an optimum. If you are interested in specific aspects of the optimization process let us know.<br>
<br>
<b>Q</b>: Did the high speed of orthAgogue result in a lower quality of the end result?<br>
<br>
<b>A</b>: No, actually the opposite is true.  Through a careful analysis of <a href='http://en.wikipedia.org/wiki/BLAST'>BLAST</a> outputs and <a href='http://www.orthomcl.org/cgi-bin/OrthoMclWeb.cgi'>OrthoMCL</a> code we've found a couple of ways to increase the resolution, which resulted in a more reliable end result. Therefore orthAgogue manages to combine high speed analysis with increased quality.<br>
<br>
<b>Q</b>: What is the difference between OrthoMCL and orthAgogue?<br>
<br>
<b>A</b>: The differences are found in (a) the programs dependencies, (b) simplicity of use,<br>
(c) sensitivity and (d) level of compatibility with the <a href='http://micans.org/mcl/index.html.bu'>MCL software</a>:<br>
<table><thead><th> <b>Aspect</b> </th><th> <b>OrthoMCL</b> </th><th> <b>orthAgogue</b> </th><th> <b>Reason</b> </th></thead><tbody>
<tr><td> Treatment of HSPs </td><td> Best only       </td><td> All               </td><td> Higher precision </td></tr>
<tr><td> Number of steps </td><td> 4               </td><td> 1                 </td><td> Ease of use   </td></tr>
<tr><td>  Need for RDBMS </td><td> Yes             </td><td> No                </td><td> Ease of use   </td></tr>
<tr><td>  FASTA header format                 </td><td> Rigid           </td><td> Flexible          </td><td> Ease of use   </td></tr>
<tr><td>  Use scores instead of e-values </td><td> No              </td><td> Yes               </td><td> Problem of '0.0' e-values </td></tr>
<tr><td>  Symmetric matrix output </td><td> No              </td><td> Yes               </td><td> Compatibility with MCL </td></tr>
<tr><td> Native MCL matrix output </td><td> No              </td><td> Yes               </td><td> Compatibility with MCL </td></tr></tbody></table>


<b>Q</b>: Why is my installation failing?<br>
<br>
<b>A</b>: Probably because orthAgogue's dependencies are not completely resolved. The dependencies are listed at the top of the installation guide <a href='https://code.google.com/p/orthagogue/wiki/InstallationGuide'>https://code.google.com/p/orthagogue/wiki/InstallationGuide</a>.<br>
<br>
<b>Q</b>: Why is my installation failing, even when the software dependencies are resolved?<br>
<br>
<b>A</b>: So far we have not experienced this. Please update us with the system description (e.g. Linux-Ubuntu, Linux-Redhat, Linux-Centos, etc.), we will then try to resolve your issue.<br>
<br>
<b>Q</b>: How do I get help on-line?<br>
<br>
<b>A</b>: Launch orthAgogue without any arguments. You'll get a synopsis of usage, a list of options explained and a number of examples of use. (We've included a softlink to the manual <a href='https://code.google.com/p/orthagogue/wiki/manual'>here</a>.)<br>
<br>
<b>Q</b>: What input file do I need?<br>
<br>
<b>A</b>: All-against-all BLAST tabular output (-m 8 option). There are almost no restrictions on the format of the <a href='http://en.wikipedia.org/wiki/FASTA_format'>FASTA file headers</a>. The only requirement is that the headers include taxon identifiers and protein identifiers unique within the taxon placed in separate fields. See command line help on the -s option.<br>
<br>
<b>Q</b>: How do I best produce my input files?<br>
<br>
<b>A</b>: Both <a href='http://blast.ncbi.nlm.nih.gov/Blast.cgi'>NCBI BLAST</a> and <a href='http://www.mpiblast.org/'>mpiBLAST</a> could be used, albeit the outputs may slightly differ.<br>
<br>
<a href='Hidden comment: 
For details about syntax highlighting, see
http://code.google.com/p/support/wiki/WikiSyntax#Code
'></a><br>
<br>
<i>An example of  mpiBLAST command:</i>
<pre><code><br>
PROJECT=cco-2011-11-17<br>
BLASTOUT=/work/mironov/outputs/mpiblast/$PROJECT<br>
FASTAFILE=/work/mironov/mpiblast/$PROJECT.fasta<br>
BLASTDB=$PROJECT.fasta # the path is specified in .ncbirc<br>
BLASTBIN=/share/apps/modulessoftware/mpiblast/mpiblast-1.6.0/bin/mpiblast<br>
mkdir $BLASTOUT<br>
mpiexec -machinefile $PBS_NODEFILE -np 96 $BLASTBIN -p blastp --debug=$BLASTOUT/$PROJECT.debug --time-profile=$BLASTOUT/$PROJECT.timeprofile --copy-via=none  -i $FASTAFILE -d $BLASTDB -e 1e-05 -o $BLASTOUT/$PROJECT.blast -m 8<br>
</code></pre>

<i>An example of NCBI BLAST command:</i>
<pre><code><br>
FASTAFILE=./$1.fasta<br>
BLASTDB=$1.fasta<br>
BLASTFILE=./$1.blast<br>
#BLASTMAT=./<br>
nohup blastall -p blastp -i $FASTAFILE -d $BLASTDB -o $BLASTFILE -e 1e-5 -m 8 -a 23  -v 100000 -b 100000 &amp;<br>
</code></pre>

<b>Q</b>: Do you have any working examples to start with?<br>
<br>
<b>A</b>: Yes, orthAgogue can be tested by issuing from the source directory:<br>
<pre><code><br>
./orthAgogue -i empirical_tests/goodProteins.blast<br>
</code></pre>

<b>Q</b>: Are there any special steps I shall take before clustering huge BlastP files?<br>
<br>
<b>A</b> Yes.<br>
<ul><li>First, remember to set the maximum memory consumption allowed by orthAgogue, i.e. using the "ulimit -Sv" command.<br>
</li><li>Second, if a single protein (in your input-blast-file) may be related to more than 50,000 other proteins (in a given taxon), then follow the instructions of the below question.</li></ul>

<b>Q</b>: What parameters shall I set when clustering huge BlastP files (e.g. a 48GB file on virtual Linux Machine with 10GB of memory)?<br>
<br>
<b>A</b>: As your file is to big to keep in memory, and orthAgogue expects having all relations for a given taxa-taxa-protein in memory, the "--disk_buffer_size" should be appropriately set. For this task, run the following script:<br>
<pre><code><br>
perl aux/estimate_disk_buffer_parameter.pl blastp_file &lt;taxon-seperator&gt; &lt;protein-index&gt;<br>
</code></pre>
After the execution (of the script), two files are created:<br>
<ul><li>at the tail of "protein_facts.txt" the "--disk_buffer_size" value is found;<br>
</li><li>in "taxa_facts.txt" we include the protein-distribution of each taxa-taxa combination.<br>
As you know have the "--disk-buffer-size" parameter, include the number (from the "protein_facts" file-tail), and then launch orthAgogue.</li></ul>

<b>Q</b>: What are the output files?<br>
<br>
<b>A</b>: By default orthAgogue generates 9 files: <i>orthologs.abc, co_orthologs.abc, inparalogs.abs, all.abc, orthologs.mcl, co_orthologs.mcl, inparalogs.mcl, all.mcl,  proteins.map</i>.<br>
<ul><li>The files with the extensions 'abc' and 'mcl' represent two different serializations of the same adjacency matrices (see the <a href='http://en.wikipedia.org/wiki/Inparanoid'>Inparanoid algorithm</a> for the definitions  of the different types of <a href='http://en.wikipedia.org/wiki/Homology_%28biology%29'>homologs</a> computed by orthAgogue).<br>
</li><li>Normally you would use either all.abc or all.mcl (all types of homologs) for clustering with MCL.<br>
</li><li>The preferred format for MCL is 'mcl', OrthoMCL uses 'abc'. For details see the <a href='http://micans.org/mcl/index.html.bu'>MCL documentation</a>.<br>
</li><li>We observed some differences in the output of MCL depending on the input file format.<br>
</li><li>The file “proteins.map” contains a mapping between the protein IDs in the BLAST file  and the protein indices used in the <code>*</code>.mcl files.</li></ul>

<b>Q</b>: Can orthAgogue emulate OrthoMCL?<br>
<br>
<b>A</b>: Yes, to emulate OrthoMCL issue from the source directory:<br>
<pre><code><br>
./orthAgogue -i empirical_tests/goodProteins.blast -b;<br>
</code></pre>
The above set of parameters is integrated in our black-box test. The starting-point of our black-box test is found in <a href='https://code.google.com/p/orthagogue/source/browse/assert_corectness.bash'>the root of the source</a>. To launch the comparison, run:<br>
<pre><code><br>
./assert_corectness.bash<br>
</code></pre>