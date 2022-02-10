__author__ = 'Sergey Osipov <Serega.Osipov@gmail.com>'


def is_child_domain_covered_by_parent_domain(parent_domain, child_domain):
    '''
    This check does not account for the distortions
    :param parent_domain:
    :param child_domain:
    :return:
    '''

    result = (min(child_domain.boundary_lons) > min(parent_domain.lon)) | (
            max(child_domain.boundary_lons) < max(parent_domain.lon)) | (
                     min(child_domain.boundary_lats) > min(parent_domain.lat)) | (
                     max(child_domain.boundary_lats) < max(parent_domain.lat))

    return result