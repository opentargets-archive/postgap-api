const resolveGeneLocation = (_, { geneId }, { geneLocationsCache }) => {
    if (geneId in geneLocationsCache) {
        const geneLocation = geneLocationsCache[geneId];
        return {
            ...geneLocation,
            chromosome: geneLocation.chr
        }
    } else {
        return null;
    }
};

export default resolveGeneLocation;
