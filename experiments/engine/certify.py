"""Additional actual-minima certificates; LLL measurements remain separate."""
import argparse
import json
import math
from pathlib import Path
import time
from experiment import (HERE, GRID, VARIANTS, build, dump, matrix, ZZ, vector,
                        norm2, name, original_constraints, member)
from shared import exact_minima as exact_engine
from coefficient_certificate import restrict


def proved_lll_xd(basis, delta):
    # The artifact's unspecified float type overflows on large ACD rows in
    # Sage 10.9. XD changes the exponent range, not the proved LLL contract.
    from fpylll.util import ReductionError
    try:
        return basis.LLL(delta=delta, algorithm='fpLLL:proved', fp='xd', transformation=True)
    except ReductionError:
        bits = max(512, 2*max(ZZ(z).nbits() for z in basis.list())+128)
        print('MPFR_LLL_FALLBACK', basis.nrows(), bits, flush=True)
        return basis.LLL(delta=delta, algorithm='fpLLL:proved', fp='rr', precision=bits, transformation=True)


exact_engine._proved_lll = proved_lll_xd
EXACT_SEARCH_LABEL = 'C++ GMP integer Bareiss exhaustive enumeration'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', type=int, required=True)
    parser.add_argument('--seconds', type=float, default=20)
    parser.add_argument('--wait', action='store_true')
    parser.add_argument('--block-size', type=int, default=0)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--recheck', action='store_true')
    parser.add_argument('--bkz-strategies',action='store_true')
    parser.add_argument('--quotient-bkz',action='store_true')
    parser.add_argument('--variant',choices=list(VARIANTS))
    parser.add_argument('--chunk-seconds',type=float,default=120)
    parser.add_argument('--prepare-only',action='store_true')
    args = parser.parse_args()
    import os
    if args.chunk_seconds<=0: raise ValueError('chunk seconds must be positive')
    os.environ['PARTITION_CHUNK_SECONDS']=str(args.chunk_seconds)
    if args.quotient_bkz:
        from bkz_precondition import precondition
        def quotient_policy(basis,delta):
            reduced,U=proved_lll_xd(basis,delta)
            if delta>0.9 and basis.nrows()>=40:
                print('QUOTIENT_BKZ_START',basis.nrows(),args.block_size,flush=True)
                reduced,V=precondition(reduced,args.block_size,use_strategies=args.bkz_strategies)
                reduced,W=proved_lll_xd(reduced,delta)
                U=W*V*U
                assert U*basis==reduced and abs(U.det())==1
                print('QUOTIENT_BKZ_DONE',basis.nrows(),flush=True)
            return reduced,U
        exact_engine._proved_lll=quotient_policy
    cfg = GRID[args.index]
    original_folder = HERE/'results'/name(cfg)
    flag_folder = HERE/'flag_results'/name(cfg)
    wrapper_folder = HERE/'wrapper_results'/name(cfg)
    flatter_folder = HERE/'flatter_results'/name(cfg)
    def ready(folder):
        return all((folder/f'{label}.json').exists() for label in VARIANTS)
    if args.wait:
        print('WAITING_FOR_MEASUREMENTS', name(cfg), flush=True)
        while not ready(original_folder) and not ready(flag_folder) and not ready(wrapper_folder) and not ready(flatter_folder):
            time.sleep(2)
    use_flag = cfg[0] == 'hm' and ready(flag_folder)
    use_flatter = not use_flag and ready(flatter_folder)
    use_wrapper = not use_flag and not use_flatter and ready(wrapper_folder)
    measurement_folder = flag_folder if use_flag else (wrapper_folder if use_wrapper else flatter_folder if use_flatter else original_folder)
    import fcntl
    original_folder.mkdir(parents=True, exist_ok=True)
    case_lock = (original_folder/'.certify.lock').open('a+')
    fcntl.flock(case_lock, fcntl.LOCK_SH if args.variant else fcntl.LOCK_EX)
    if args.variant:
        variant_lock=(original_folder/('.certify-'+args.variant+'.lock')).open('a+')
        fcntl.flock(variant_lock,fcntl.LOCK_EX)
    data = build(cfg)
    constraints = original_constraints(data)
    for label in VARIANTS:
        if args.variant and label!=args.variant: continue
        source = measurement_folder / f'{label}.json'
        if not source.exists():
            raise FileNotFoundError(source)
        doc = json.loads(source.read_text())
        target = original_folder / (label + '.certificate.json')
        if target.exists() and json.loads(target.read_text()).get('exact') and not args.recheck:
            continue
        basis = matrix(ZZ, doc['reduced_rows'])
        scale = ZZ(doc['common_coordinate_divisor'])
        reduced_scale = basis.apply_map(lambda z: z // scale)
        start = time.monotonic()
        if doc['exact_by_divisibility']:
            if args.prepare_only:
                print('NO_PREPARATION_NEEDED',name(cfg),label,flush=True)
                continue
            witnesses = [list(map(str, basis.row(i))) for i in doc['shortest_row_indices']]
            certificate = dict(exact=True, method='coordinate_divisibility_lower_bound',
                lambda1_sq=doc['lambda1_sq'], lambda2_sq=doc['lambda2_sq'],
                lambda1_witness=witnesses[0], lambda2_witness=witnesses[1])
        else:
            shortest_index = doc['shortest_row_indices'][0]
            first_sq = norm2(reduced_scale.row(shortest_index))
            first_basis, first_exclusion = restrict(reduced_scale, first_sq)
            first_exact = first_basis.nrows() == 1
            upper_sq = max(norm2(reduced_scale.row(i)) for i in doc['shortest_row_indices'])
            restricted, exclusion = restrict(reduced_scale, upper_sq)
            assert restricted.nrows() >= 2
            bkz_tag=(('-qbkz' if args.quotient_bkz else '-bkz')+str(args.block_size)+('s' if args.bkz_strategies else '')) if args.block_size else ''
            checkpoint = HERE / 'certificates' / name(cfg) / (label + ('-flag' if use_flag else '-wrapper' if use_wrapper else '-flatter' if use_flatter else '') + bkz_tag)
            dump(checkpoint / 'coefficient_exclusion.json', exclusion)
            dump(checkpoint / 'lambda1_exclusion.json', first_exclusion)
            print('EXCLUSION', name(cfg), label, basis.nrows(), '->', restricted.nrows(), flush=True)
            if args.block_size and restricted.nrows() >= 40 and not args.quotient_bkz:
                from bkz_precondition import precondition
                before = time.monotonic()
                preconditioning = checkpoint / 'preconditioning.json'
                if preconditioning.exists():
                    saved = json.loads(preconditioning.read_text())
                    assert saved['block_size'] == args.block_size
                    transform = matrix(ZZ, saved['transform_rows'])
                    assert abs(transform.det()) == 1
                    restricted = transform*restricted
                    print('BKZ_REUSED', name(cfg), label, flush=True)
                else:
                    print('BKZ_START', name(cfg), label, args.block_size, flush=True)
                    restricted, transform = precondition(restricted, args.block_size,use_strategies=args.bkz_strategies)
                    dump(preconditioning, dict(block_size=args.block_size,
                         loops=4, precision=512, transform_rows=[list(map(str,row)) for row in transform.rows()],
                         seconds=time.monotonic()-before))
                    print('BKZ_DONE', name(cfg), label, time.monotonic()-before, flush=True)
            if args.prepare_only:
                exact_engine.prepare_successive_minima(restricted,checkpoint,workers=1)
                print('PREPARED_ONLY',name(cfg),label,str(checkpoint),flush=True)
                continue
            engine = exact_engine.exact_successive_minima(restricted,
                checkpoint,
                workers=args.workers, memory_gb=4, max_seconds=args.seconds,
                chunk_nodes=20000, chunk_seconds=1, progress_seconds=10)
            engine['actual_enumeration_backend'] = EXACT_SEARCH_LABEL
            engine['quotient_bkz']=dict(block_size=args.block_size,strategies=args.bkz_strategies) if args.quotient_bkz else None
            engine['enumeration_backend'] = EXACT_SEARCH_LABEL
            dump(checkpoint/'result.json',engine)
            dump(checkpoint/'runner_result.json',engine)
            certificate = dict(exact=engine['exact'], method='exhaustive_integer_enumeration',
                status=engine['status'], lambda1_sq=None, lambda2_sq=None,
                engine_checkpoint=engine['checkpoint'],
                ambient_rank=basis.nrows(), certified_search_rank=restricted.nrows(),
                exclusion_upper_sq=str(upper_sq*scale**2))
            certificate['enumeration_backend'] = EXACT_SEARCH_LABEL
            if first_exact and engine.get('lambda1_sq') is None:
                certificate.update(lambda1_sq=str(first_sq*scale**2),
                    lambda1_witness=list(map(str,basis.row(shortest_index))),
                    lambda1_method='one_remaining_coefficient_direction')
            for i in (1, 2):
                if engine.get(f'lambda{i}_sq') is not None:
                    certificate[f'lambda{i}_sq'] = str(ZZ(engine[f'lambda{i}_sq']) * scale**2)
                    certificate[f'lambda{i}_witness'] = [str(ZZ(z) * scale) for z in engine[f'lambda{i}_witness']]
        for i in (1, 2):
            if certificate.get(f'lambda{i}_sq') is not None:
                witness = vector(ZZ, certificate[f'lambda{i}_witness'])
                assert norm2(witness) == ZZ(certificate[f'lambda{i}_sq'])
                assert member(witness, data, constraints, VARIANTS[label])
                certificate[f'log2_lambda{i}'] = math.log2(int(certificate[f'lambda{i}_sq']))/2
        if certificate['exact']:
            assert matrix(ZZ, [certificate['lambda1_witness'], certificate['lambda2_witness']]).rank() == 2
            assert ZZ(certificate['lambda1_sq']) <= ZZ(certificate['lambda2_sq'])
            if certificate['method'] == 'exhaustive_integer_enumeration':
                assert ZZ(certificate['lambda2_sq']) <= upper_sq*scale**2
        certificate.update(case=name(cfg), variant=label, common_scale=str(scale),
                           seconds=time.monotonic()-start, lll_float_type='xd',
                           lambda1_exact=certificate.get('lambda1_sq') is not None,
                           lambda2_exact=bool(certificate['exact']))
        certificate['measurement_file'] = str(source)
        dump(target, certificate)
        print('CERTIFIED', name(cfg), label, certificate['exact'],
              certificate.get('log2_lambda1'), certificate.get('log2_lambda2'), flush=True)


if __name__ == '__main__':
    main()
