import axios from 'axios';

const ENSEMBL_API_BASE = 'https://rest.ensembl.org/';
const ENSEMBL_API_VARIATION = 'variation/homo_sapiens';
const ENSEMBL_API_LOOKUP = 'lookup/id';

function transformEnsemblGene(d) {
    const {
      id,
      display_name,
      description,
      start,
      end,
      strand,
      seq_region_name,
      biotype,
      Transcript
    } = d;
    let canonicalTranscript = Transcript.filter(t => t.is_canonical === 1).map(
      t => {
        const { id, start, end, strand, Exon, Translation } = t;
        const exons = Exon.map(ex => ({
          id: ex.id,
          start: ex.start,
          end: ex.end
        }));
        const translation = Translation
          ? {
              translationStart: Translation.start,
              translationEnd: Translation.end
            }
          : {};
        const tss = strand === 1 ? start : end; // tss depends on strand
        return {
          id,
          start,
          end,
          forwardStrand: (strand === 1),
          exons,
          tss,
          ...translation
        };
      }
    );
    if (canonicalTranscript.length === 1) {
      canonicalTranscript = canonicalTranscript[0];
    } else {
      canonicalTranscript = null; // no transcript
    }
    return {
      id,
      description,
      start,
      end,
      forwardStrand: (strand === 1),
      biotype,
      name: display_name,
      chromosome: seq_region_name,
      canonicalTranscript
    };
  }
  
  function transformEnsemblVariant(d) {
    return {
      id: d.name,
      maf: d.MAF,
      ancestralAllele: d.ancestral_allele,
      minorAllele: d.minor_allele,
      position: d.mappings[0].start,
      chromosome: d.mappings[0].seq_region_name
    };
  }
  

export default {
    fetchVariants(variantIds) {
      const url = `${ENSEMBL_API_BASE}${ENSEMBL_API_VARIATION}`;
      const body = {
        ids: variantIds,
      };
      return axios.post(url, body).then(response => {
        console.log(`ensemblClient.fetchVariants: ${response.status} (${variantIds.length} variants)`);
        const variantsRaw = response.data;
        const variantsClean = {};
        for (const [variantId, variantObj] of Object.entries(variantsRaw)) {
            variantsClean[variantId] = transformEnsemblVariant(variantObj);
        }
        return variantsClean;
      }).catch(function (error) {
        console.log(error);
      });
    },
    fetchGenes(geneIds) {
      const url = `${ENSEMBL_API_BASE}${ENSEMBL_API_LOOKUP}`;
      const body = {
        ids: geneIds,
        expand: true,
      };
      return axios.post(url, body).then(response => {
        console.log(`ensemblClient.fetchGenes: ${response.status} (${geneIds.length} genes)`);
        const genesRaw = response.data;
        const genesClean = {};
        for (const [geneId, geneObj] of Object.entries(genesRaw)) {
            genesClean[geneId] = transformEnsemblGene(geneObj);
        }
        return genesClean;
      }).catch(function (error) {
        console.log(error);
      });
    },
    // fetchSearch(query) {
    //   const url = `${ENSEMBL_API_BASE}${ENSEMBL_API_VARIATION}`;
    //   const body = {
    //     ids: [query],
    //   };
    //   return axios.post(url, body).then(response => {
    //     const data = [];
    //     if (response.data && response.data[query]) {
    //       data.push({
    //         id: query,
    //         name: query,
    //         type: 'variant',
    //       });
    //     }
    //     return data;
    //   });
    // },
  };