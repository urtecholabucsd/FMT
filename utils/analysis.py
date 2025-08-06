import pandas as pd
import os
import warnings
from micom.qiime_formats import load_qiime_medium
from collections import Counter
from micom import interaction
from micom.workflows import build, save_results, load_results, GrowthResults, grow
import numpy as np
from sklearn.metrics.pairwise import cosine_distances, euclidean_distances
from typing import List, Any, Union, Dict, Tuple, Optional, Literal
from sklearn.metrics.pairwise import pairwise_distances
from scipy.stats import mannwhitneyu
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings('ignore')

def custom_sort(taxa: List[str]) -> List[str]:
    """
    Custom sort function for sorting the taxa

    Acceptable formats: 

    - String Format: k__Bacteria, p__Firmicutes, c__Clostridia, o__Clostridiales, f__Lachnospiraceae, g__Blautia
    - String and Numeric Format: AZC3J_1, AZC3J_2, AZC3J_3, AZC3J_4, AZC3J_5, AZC3J_6

    If both formats are detected in the input list, raises a ValueError.
    If only string format is provided, the list is sorted lexicographically.
    If string and numeric format is provided, it sorts by the prefix (before '_') 
    and then by the numeric part.

    :param taxa: list of taxa
    """
    string_only = all('_' not in x for x in taxa)
    string_and_numeric = all('_' in x and x.split('_')[1].isdigit() for x in taxa)

    if not (string_only or string_and_numeric):
        raise ValueError('Taxa list contains mixed types or invalid labels')
    if string_only:
        return sorted(taxa)
    elif string_and_numeric:
        return sorted(taxa, key=lambda x: (x.split('_')[0], int(x.split('_')[1])))


def _mes(df: pd.DataFrame) -> pd.Series:
    """Helper to calculate the MES score."""
    cn = Counter(df.direction)
    p, c = cn["export"], cn["import"]
    return pd.Series(2.0 * p * c / (p + c), index=["MES"])


def _log2foldchange(df: pd.DataFrame):
    d0_abundance = df[df['Day'] == 'D0']['abundance'].values[0]
    d5_abundance = df[df['Day'] == 'D5']['abundance'].values[0]
    return np.log2((d5_abundance + 1) / (d0_abundance + 1))


def harmonize_matrices(mat1, mat2):
    # Unify metabolites (columns)
    combined_cols = mat1.columns.union(mat2.columns)
    mat1 = mat1.reindex(columns=combined_cols, fill_value=0)
    mat2 = mat2.reindex(columns=combined_cols, fill_value=0)

    # Unify taxa (rows)
    combined_rows = mat1.index.union(mat2.index)
    mat1 = mat1.reindex(index=combined_rows, fill_value=0)
    mat2 = mat2.reindex(index=combined_rows, fill_value=0)

    return mat1, mat2


def calculate_metabolic_interactions(
        res: GrowthResults,
        out_folder: str,
        taxa: List[str] = None,
        threads: int = 11) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Calculate metabolic interactions and metabolic exchange score between the microbes in the dataset
    :param res:
    :param out_folder:
    :param metric:
    :param taxa:
    :param threads:
    :return:
    """

    if 'interaction_results.csv' in os.listdir(out_folder):
        interaction_res = pd.read_csv(os.path.join(out_folder, 'interaction_results.csv'))
        print('Interaction results already exist, loading the results')
    else:
        interaction_res = interaction.interactions(res, taxa=taxa, threads=threads)
        interaction_res.to_csv(os.path.join(out_folder, 'interaction_results.csv'), index=False)

    cutoff = res.exchanges.tolerance
    fluxes = res.exchanges[
        (res.exchanges.flux.abs() > cutoff) &
        (res.exchanges.taxon != 'medium')
        ]
    if 'mes.csv' in os.listdir(out_folder):
        mes = pd.read_csv(os.path.join(out_folder, 'mes.csv'))
        print('MES scores already exist, loading the results')
    else:
        mes = fluxes.groupby(["metabolite", "sample_id"]).apply(_mes).reset_index()
        mes = mes.merge(res.annotations.drop_duplicates(subset=['metabolite']),
                        on='metabolite', how='inner')
        mes.to_csv(os.path.join(out_folder, 'mes.csv'), index=False)

    return mes, interaction_res


def feature_extraction(
        distance_df_between: Optional[pd.DataFrame] = None,
        distance_df_within: Optional[pd.DataFrame] = None,
        enrichment: Optional[pd.DataFrame] = None,
        N_values: List[int] = [5, 10, 15, 20]
) -> pd.DataFrame:
    """

    :param distance_df_between:
    :param distance_df_within:
    :param enrichment:
    :param N_values:
    :return:
    """
    if distance_df_between is None and distance_df_within is None:
        raise ValueError('Please provide the distance data for calculating the feature importance')

    if enrichment is not None:
        if 'colonization' in enrichment.columns:
            if distance_df_between is not None:
                distance_df_between = distance_df_between.reset_index().merge(enrichment, on='taxon',
                                                                              how='left').set_index('taxon')
            if distance_df_within is not None:
                distance_df_within = distance_df_within.reset_index().merge(enrichment, on='taxon',
                                                                            how='left').set_index('taxon')
        else:
            raise ValueError('Please provide the enrichment data with colonization outcomes')
    else:
        if distance_df_between is not None and 'colonization' not in distance_df_between.columns:
            raise ValueError('Please provide the enrichment data with colonization outcomes')
        if distance_df_within is not None and 'colonization' not in distance_df_within.columns:
            raise ValueError('Please provide the enrichment data with colonization outcomes')

    final = []
    final_results = []
    for N in N_values:
        if distance_df_between is not None:
            average_nearest_N_distances = []
            for row in distance_df_between.values:
                nearest_N = np.sort(row[:-1])[:N]
                average_distance = np.mean(nearest_N)
                average_nearest_N_distances.append(average_distance)

            results_df = pd.DataFrame({
                'taxon': distance_df_between.index,
                'average_distance': average_nearest_N_distances,
                'colonization': distance_df_between['colonization'],
                'N': N,
                'comparisons': 'between'
            })
            final_results.append(results_df)

        if distance_df_within is not None:
            average_nearest_N_distances = []
            for row in distance_df_within.values:
                nearest_N = np.sort(row[:-1])[1:N + 1]  # avoid self-interactions
                average_distance = np.mean(nearest_N)
                average_nearest_N_distances.append(average_distance)

            results_df = pd.DataFrame({
                'taxon': distance_df_within.index,
                'average_distance': average_nearest_N_distances,
                'colonization': distance_df_within['colonization'],
                'N': N,
                'comparisons': 'within'
            })
            final_results.append(results_df)

        final = pd.concat(final_results, ignore_index=True)

    for N in N_values:
        for comp in final['comparisons'].unique():
            subset = final[(final['N'] == N) & (final['comparisons'] == comp)]
            colonizers = subset[subset['colonization'] == 'colonized']
            non_colonizers = subset[subset['colonization'] == 'non-colonized']
            if len(non_colonizers) == 0:
                print(f"No non-colonized MAGs for N={N} and comparison={comp}")
                continue
            elif len(colonizers) == 0:
                print(f"No colonized MAGs for N={N} and comparison={comp}")
                continue
            _, p_val = mannwhitneyu(colonizers['average_distance'], non_colonizers['average_distance'])

            print(f"For N={N}:")
            print(f"Mean average distance for colonized MAGs: {colonizers['average_distance'].mean()}")
            print(f"Mean average distance for non-colonized MAGs: {non_colonizers['average_distance'].mean()}")
            print(f"p-value: {p_val}\n")

    return final


def metabolic_exchange_distance(exchanges_df_t1: pd.DataFrame,
                                donor: str,
                                recipient: str,
                                exchanges_df_t2: pd.DataFrame,
                                metabolites: Optional[List[str]] = None,
                                metric: Literal['euclidean', 'cosine', 'mahalanobis'] = 'euclidean',
                                direction: Literal['import', 'export', 'delta'] = 'import') -> pd.DataFrame:
    """

    :param exchanges_df_t2:
    :param groups:
    :param exchanges_df_t1:
    :param metric:
    :param direction:
    :return:
    """

    if direction not in ['import', 'export', 'delta']:
        raise ValueError('Please provide a valid direction for the metabolic exchange')

    if metabolites is None:
        metabolites = pd.concat([exchanges_df_t1['metabolite'], exchanges_df_t2['metabolite']], axis=0).unique()

    exchanges_t1 = exchanges_df_t1[['taxon', 'sample_id', 'flux', 'metabolite', 'direction']].copy()
    exchanges_t1 = exchanges_t1[exchanges_t1['taxon'].str.contains(donor)]

    exchanges_t2 = exchanges_df_t2[['taxon', 'sample_id', 'flux', 'metabolite', 'direction']].copy()
    exchanges_t2 = exchanges_t2[exchanges_t2['taxon'].str.contains(recipient)]

    taxa = pd.concat([exchanges_t1['taxon'], exchanges_t2['taxon']], axis=0).unique()
    taxa = custom_sort(taxa)
    if direction == 'delta':
        def _apply_delta(df):
            temp_import = df[df['direction'] == 'import']
            temp_export = df[df['direction'] == 'export']

            temp_import = temp_import.pivot_table(index='taxon',
                                                  columns='metabolite',
                                                  values='flux',
                                                  fill_value=0.0,
                                                  aggfunc='mean')

            temp_export = temp_export.pivot_table(index='taxon',
                                                  columns='metabolite',
                                                  values='flux',
                                                  fill_value=0.0,
                                                  aggfunc='mean')

            temp_import = temp_import.reindex(index=taxa, columns=metabolites, fill_value=0.0)
            temp_export = temp_export.reindex(index=taxa, columns=metabolites, fill_value=0.0)

            temp = temp_export - temp_import
            return temp

        exchanges_t1 = _apply_delta(exchanges_t1)
        exchanges_t2 = _apply_delta(exchanges_t2)

    else:
        exchanges_t1 = exchanges_t1[exchanges_t1['direction'] == direction]
        exchanges_t1 = exchanges_t1.pivot_table(index='taxon',
                                                columns='metabolite',
                                                values='flux',
                                                fill_value=0.0,
                                                aggfunc='mean')

        exchanges_t1 = exchanges_t1.reindex(index=taxa, columns=metabolites, fill_value=0.0)

        exchanges_t2 = exchanges_t2[exchanges_t2['direction'] == direction]
        exchanges_t2 = exchanges_t2.pivot_table(index='taxon',
                                                columns='metabolite',
                                                values='flux',
                                                fill_value=0.0,
                                                aggfunc='mean')

        exchanges_t2 = exchanges_t2.reindex(index=taxa, columns=metabolites, fill_value=0.0)

    if metric == 'euclidean':
        distance_matrix = euclidean_distances(exchanges_t1.values, exchanges_t2.values)
    elif metric == 'cosine':
        distance_matrix = cosine_distances(exchanges_t1.values, exchanges_t2.values)
    elif metric == 'mahalanobis':
        combined = pd.concat([exchanges_t1, exchanges_t2], axis=0)
        from sklearn.covariance import LedoitWolf
        cov = LedoitWolf().fit(combined)
        S_inv = cov.precision_
        distance_matrix = pairwise_distances(exchanges_t1.values, exchanges_t2.values, metric='mahalanobis', VI=S_inv)
    else:
        raise ValueError('Please provide a valid metric for calculating the distance')

    distance_df = pd.DataFrame(distance_matrix, index=exchanges_t1.index, columns=exchanges_t2.index)
    return distance_df


def metabolic_interaction_distance(donor_df: pd.DataFrame,
                                   donor: str,
                                   recipient: str,
                                   recipient_df: pd.DataFrame,
                                   metabolites: Optional[List[str]] = None,
                                   method: Literal['euclidean', 'cosine', 'mahalanobis'] = 'euclidean') -> pd.DataFrame:

    if method not in ['euclidean', 'cosine', 'mahalanobis']:
        raise ValueError('Please provide a valid method for calculating the distance')

    if not all([col in donor_df.columns \
                for col in ['metabolite', 'focal', 'partner', 'sample_id', 'flux', 'class']]) and \
            not all([col in recipient_df.columns \
                     for col in ['metabolite', 'focal', 'partner', 'sample_id', 'flux', 'class']]):
        raise ValueError('Please provide a valid interaction dataframe with metabolite, focal, partner, sample_id, '
                         'flux, and class columns')

    int_donor = donor_df.groupby(['metabolite',
                                        'focal', 'partner', 'class'])['flux'].mean().to_frame().reset_index().rename(
        columns={'focal': 'taxon'})
    int_recipient = recipient_df.groupby(['metabolite',
                                        'focal', 'partner', 'class'])['flux'].mean().to_frame().reset_index().rename(
        columns={'focal': 'taxon'})

    if metabolites is None:
        metabolites = sorted(pd.concat([int_donor['metabolite'], int_recipient['metabolite']], axis=0).unique())


    int_donor = int_donor.pivot_table(index='taxon',
                                      columns='metabolite',
                                      values='flux',
                                      fill_value=0.0,
                                      aggfunc='sum')
    
    int_donor = int_donor.reindex(columns=metabolites, fill_value=0.0)

    int_recipient = int_recipient.pivot_table(index='taxon',
                                              columns='metabolite',
                                              values='flux',
                                              fill_value=0.0,
                                              aggfunc='sum')
    
    int_recipient = int_recipient.reindex(columns=metabolites, fill_value=0.0)


    if method == 'euclidean':
        distance_matrix = euclidean_distances(int_donor.values, int_recipient.values)
    elif method == 'cosine':
        distance_matrix = cosine_distances(int_donor.values, int_recipient.values)
    elif method == 'mahalanobis':
        combined = pd.concat([int_donor, int_recipient], axis=0)
        from sklearn.covariance import LedoitWolf
        cov = LedoitWolf().fit(combined)
        S_inv = cov.precision_
        distance_matrix = pairwise_distances(int_donor.values,
                                             int_recipient.values,
                                             metric='mahalanobis', VI=S_inv)
    else:
        raise ValueError('Please provide a valid method for calculating the distance')

    dist_df = pd.DataFrame(distance_matrix, index=int_donor.index, columns=int_recipient.index)
    return dist_df


def sample_level_distance(combined_exchange_dfs: Dict[str, pd.DataFrame],
                          method: Literal['euclidean', 'cosine', 'jaccard'] = 'euclidean',
                          direction: Literal['import', 'export', 'delta'] = 'delta') -> pd.DataFrame:

    if method not in ['euclidean', 'cosine', 'jaccard']:
        raise ValueError('Please provide a valid method for calculating the distance')

    if direction not in ['import', 'export', 'delta']:
        raise ValueError('Please provide a valid direction for the metabolic exchange')

    metabolites = set()
    for df in combined_exchange_dfs.values():
        metabolites.update(df['metabolite'].unique())

    if direction == 'delta':
        ex = {}
        for key, df in combined_exchange_dfs.items():
            mags = sorted(df['taxon'].unique())
            temp_import = df[df['direction'] == 'import']
            temp_export = df[df['direction'] == 'export']

            temp_import = temp_import.pivot_table(index='taxon',
                                         columns='metabolite',
                                         values='flux',
                                         fill_value=0.0,
                                         aggfunc='mean')

            temp_export = temp_export.pivot_table(index='taxon',
                                         columns='metabolite',
                                         values='flux',
                                         fill_value=0.0,
                                         aggfunc='mean')

            temp_import = temp_import.reindex(index=mags, columns=metabolites, fill_value=0.0)
            temp_export = temp_export.reindex(index=mags, columns=metabolites, fill_value=0.0)

            temp = temp_export - temp_import
            temp = temp.sum(axis=0).to_frame().T.rename(index={0: key})
            ex[key] = temp

    elif direction == 'import':
        ex = {}
        for key, df in combined_exchange_dfs.items():
            mags = sorted(df['taxon'].unique())
            temp = df[df['direction'] == 'import']

            temp = temp.pivot_table(index='taxon',
                                    columns='metabolite',
                                    values='flux',
                                    fill_value=0.0,
                                    aggfunc='mean')

            temp = temp.reindex(index=mags, columns=metabolites, fill_value=0.0)
            temp = temp.sum(axis=0).to_frame().T.rename(index={0: key})
            ex[key] = temp

    elif direction == 'export':
        ex = {}
        for key, df in combined_exchange_dfs.items():
            mags = sorted(df['taxon'].unique())
            temp = df[df['direction'] == 'export']

            temp = temp.pivot_table(index='taxon',
                                    columns='metabolite',
                                    values='flux',
                                    fill_value=0.0,
                                    aggfunc='mean')

            temp = temp.reindex(index=mags, columns=metabolites, fill_value=0.0)
            temp = temp.sum(axis=0).to_frame().T.rename(index={0: key})
            ex[key] = temp

    else:
        raise ValueError('Please provide a valid direction for the metabolic exchange')

    ex = pd.concat(ex.values(), axis=0)
    ex_scaled = StandardScaler().fit_transform(ex.values)
    ex_scaled = pd.DataFrame(ex_scaled, index=ex.index, columns=ex.columns)

    distance_matrix = pdist(ex_scaled.values, metric=method)

    distance_matrix = squareform(distance_matrix)
    distance_df = pd.DataFrame(distance_matrix, index=ex.index, columns=ex.index)

    return distance_df

def metabolic_distribution(df: pd.DataFrame,
                           regularization: float = 1e-6) -> pd.DataFrame:
    """
    Compute the similarity between the distribution of each metabolite and the null distribution.
    :param df:
    :param regularization:
    :return:
    """

    res = {"metabolite": [], "mahalanobis_distance": []}

    for col in df.columns:
        col_vec = df[col].values
        null_vec = np.zeros_like(col_vec)  # Null vector (reference)

        cov = np.cov(col_vec, rowvar=False)

        cov_m = np.atleast_2d(cov) + np.eye(1) * regularization  # Regularization to avoid singular matrix

        if np.isclose(np.linalg.det(cov_m), 0):
            res['mahalanobis_distance'].append(0)
        else:
            S_inv = np.linalg.inv(cov_m)
            delta = col_vec - null_vec
            dist = np.sqrt(np.sum(delta * delta / np.diag(S_inv)))
            res['mahalanobis_distance'].append(dist)
        res['metabolite'].append(col)

    return pd.DataFrame(res)


if __name__ == '__main__':
    pass
