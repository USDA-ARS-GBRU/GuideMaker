<h3
	style="color: #003087"
	>Pooled CRISPR Experiments
</h3>
<div
	style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>Experiments that target the entire genome, or many genes at once, are typically performed in pooled experiments where 100-100,000+ targets are tested simultaneously. The pooled oligonucleotides for each gRNA are cloned in one batch and used simultaneously in the designed experiment. Each gRNA sequence acts as a barcode that can be quantified with high-throughput sequencing to elucidate each target's relative importance under the experimental conditions.
</div>
<h3 style="color: #003087" >Vectors for gRNA Cloning</h3>
<div
	style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>Genome-scale CRISPR experiments require a gRNA vector amenable to high-throughput cloning, most often through [Golden Gate cloning](https://blog.addgene.org/plasmids-101-golden-gate-cloning), a restriction enzyme-dependent reaction. Plasmids to express gRNA are available from Addgene can be found at the link below, though not all of these are compatible with high-throughput cloning.  
</div>
<h3 style="color: #003087" >Addgene: CRISPR Plasmids - Empty gRNA Vectors</h3>
<div
style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>After running GuideMaker, the designed gRNA output can be downloaded and with minor adjustments, the targets can be ordered as oligos for cloning. Pooled oligonucleotides can be purchased from several vendors, including those listed below. Pool sizes vary from 100 to over 200,000 oligonucleotides. Vendor specifications for the number of oligos, oligo length, and cost per bp vary widely. For bacterial genome-scale experiments, as of 2021, Genscript offers pool sizes of 12,472 and 91,766 with up to 79 bp per oligo for list prices of $1600 and $4,000, respectively.

Some example vendors are:

*   [GenScript](https://www.genscript.com/precise-synthetic-oligo-pools.html)
*   [Twist Bioscience](https://www.twistbioscience.com/products/oligopools)
*   [Agilent](https://www.agilent.com/en/product/sureprint-oligonucleotide-library-synthesis/oligonucleotide-library-synthesis/sureprint-oligonucleotide-libraries-288039)
*   [arbor biosciences](https://arborbiosci.com/oligos-and-arrays/dna-and-rna-oligo-pools/)
</div>
<div
style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>Most pools require amplification before cloning to convert the ssDNA to dsDNA and increase the concentration for efficient cloning. Accordingly, adding a constant region at the 3' end for primer binding is recommended. Sub-pools can also be amplified by adding unique constant regions to some oligos, enabling the large-scale synthesis to be split amongst organisms or specific targets in a single organism. Because Golden Gate cloning utilizes restriction enzymes, filtering gRNA designs with the cognate restriction enzyme recognition sites is necessary, a feature found in GuideMaker. A general protocol for cloning pooled gRNA from synthesized oligonucleotides from IDT is linked below, though similar workflows can be used for pools from other vendors.

*   [Cloning high-quality CRISPR libraries with oPools Oligo Pools (SYB-10182-PR12/2019)](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/user-submitted-method/cloning-high-quality-crispr-libraries-with-opools-oligo-pools-user-method.pdf?sfvrsn=3db31607_7)
*   [Addgene: Guide to Using Pooled Libraries](https://www.addgene.org/guides/pooled-libraries/)
 </div>
<h3 style="color: #003087">Pooled CRISPR Data Analysis</h3>
<div
style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>After the experiment, the cells are collected and DNA is isolated. The target sequence is then amplified and adaptors for high-throughput sequencing added. Several data analysis pipelines have been developed to identify target sequences over-represented or under-represented in the pool. The manuscript by Wang et al. (2019) provides a protocol for using a high-quality tool with these capabilities. 
 </div>
<h3 style="color: #003087">Citation</h3>
<div
style="font-family:Hoefler Text;
	font-size: 14px;
	color: "black"
	>Wang, B., Wang, M., Zhang, W. et al. Integrative analysis of pooled CRISPR genetic screens using MAGeCKFlute. Nat Protoc 14, 756–780 (2019). https://doi.org/10.1038/s41596-018-0113-7
 </div>