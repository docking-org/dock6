/*                                                                    */
/*                        Copyright UCSF, 1997                        */
/*                                                                    */

#include "define.h"
#include "utility.h"
#include "global.h"
#include "search.h"


/* ///////////////////////////////////////////////////////////// */

void allocate_search (SEARCH *search)
{
  int i;

  if (search->max_size > 0)
  {
    ecalloc
    (
      (void **) &search->count,
      search->max_size + 1,
      sizeof (int),
      "search count",
      global.outfile
    );

    ecalloc
    (
      (void **) &search->total,
      search->max_size + 1,
      sizeof (int),
      "search total",
      global.outfile
    );

    ecalloc
    (
      (void **) &search->target,
      search->max_size + 1,
      sizeof (int *),
      "search target",
      global.outfile
    );

    ecalloc
    (
      (void **) &search->origin,
      search->max_size + 1,
      sizeof (int *),
      "search origin",
      global.outfile
    );

    for (i = 0; i < search->max_size + 1; i++)
    {
      ecalloc
      (
        (void **) &search->target[i],
        search->max_size,
        sizeof (int),
        "search target",
        global.outfile
      );

      ecalloc
      (
        (void **) &search->origin[i],
        search->max_size,
        sizeof (int),
        "search origin",
        global.outfile
      );
    }

    ecalloc
    (
      (void **) &search->neighbor_id,
      search->max_size + 1,
      sizeof (int),
      "search neighbor id",
      global.outfile
    );

    ecalloc
    (
      (void **) &search->log,
      search->max_size + 1,
      sizeof (int),
      "search log",
      global.outfile
    );
  }
}

/* ///////////////////////////////////////////////////////////// */

void reset_search (SEARCH *search)
{
  int i;

  if (search->max_size > 0)
  {
    memset (search->count, 0, (search->max_size + 1) * sizeof (int));
    memset (search->total, 0, (search->max_size + 1) * sizeof (int));

    for (i = 0; i < search->max_size + 1; i++)
    {
      memset (search->target[i], 0, search->max_size * sizeof (int));
      memset (search->origin[i], 0, search->max_size * sizeof (int));
    }

    memset (search->neighbor_id, 0, (search->max_size + 1) * sizeof (int));
    memset (search->log, 0, (search->max_size + 1) * sizeof (int));
  }

  search->complete_flag = FALSE;
  search->radius = 0;
  search->total_size = 0;
}

/* ///////////////////////////////////////////////////////////// */

void free_search (SEARCH *search)
{
  int i;

  efree ((void **) &search->count);
  efree ((void **) &search->total);

  if (search->max_size > 0)
  {
    for (i = 0; i < search->max_size + 1; i++)
    {
      efree ((void **) &search->target[i]);
      efree ((void **) &search->origin[i]);
    }
  }

  efree ((void **) &search->target);
  efree ((void **) &search->origin);
  efree ((void **) &search->neighbor_id);
  efree ((void **) &search->log);

  search->max_size = 0;
}

/* ///////////////////////////////////////////////////////////// */

void reallocate_search (SEARCH *search)
{
  if (search->total_size > search->max_size)
  {
    free_search (search);
    search->max_size = search->total_size;
    allocate_search (search);
  }
}


/* //////////////////////////////////////////////////////////////////////

Subroutine to traverse a list of linked nodes non-recursively, breadth-first.

11/96 te

////////////////////////////////////////////////////////////////////// */

int breadth_search
(
  SEARCH	*search,
  void		*nodes,
  int		total,
  int		get_neighbor (void *, int, int),
  void		flag_neighbor (void *, int, int, int),
  int		*seed,
  int		seed_total,
  int		avoid,
  int		iteration
)
{
  int i;
  int node;
  int neighbor;

  if (iteration == 0)
  {
    reset_search (search);
    search->total_size = total;
    reallocate_search (search);

    if ((avoid >= 0) && (avoid < total))
      search->log[avoid] = TRUE;

    for (i = 0; i < seed_total; i++)
    {
      if ((seed[i] >= 0) && (seed[i] < total))
      {
        search->log[seed[i]] = TRUE;
        search->target[search->radius][search->total[search->radius]] = seed[i];
        search->origin[search->radius][search->total[search->radius]] = NEITHER;
        search->total[search->radius]++;
      }

      else
        exit (fprintf (global.outfile,
          "ERROR breadth_search: seed value inappropriate (%d).\n", seed[i]));
    }
  }

  else
    search->count[search->radius]++;

  if (search->count[search->radius] < search->total[search->radius])
    return search->target[search->radius][search->count[search->radius]];

  else
  {
    for
    (
      search->count[search->radius] = 0;
      search->count[search->radius] < search->total[search->radius];
      search->count[search->radius]++
    )
    {
      node = search->target[search->radius][search->count[search->radius]];

      for
      (
        search->neighbor_id[search->radius] = 0;
        (neighbor = get_neighbor
          (nodes, node, search->neighbor_id[search->radius])) != EOF;
        search->neighbor_id[search->radius]++
      )
      {
        if (flag_neighbor)
          flag_neighbor
          (
            nodes,
            node,
            search->neighbor_id[search->radius],
            !search->log[neighbor]
          );

        if (!search->log[neighbor])
        {
          search->target
            [search->radius + 1]
            [search->total[search->radius + 1]] = neighbor;
          search->origin
            [search->radius + 1]
            [search->total[search->radius + 1]] = node;

          search->total[search->radius + 1]++;
          search->log[neighbor] = TRUE;
        }
      }
    }

    if (search->total[search->radius + 1] > 0)
    {
      search->radius++;
      return search->target[search->radius][search->count[search->radius]];
    }

    else
    {
      search->complete_flag = TRUE;
      return EOF;
    }
  }
}

/* //////////////////////////////////////////////////////////////////////

Subroutine to traverse a list of linked nodes non-recursively, depth-first.

11/96 te

////////////////////////////////////////////////////////////////////// */

int depth_search
(
  SEARCH	*search,
  void		*nodes,
  int		total,
  int		get_neighbor (void *, int, int),
  void		flag_neighbor (void *, int, int, int),
  int		*seed,
  int		seed_total,
  int		avoid,
  int		iteration
)
{
  int i;
  int node;
  int neighbor;
  int terminus_found;

  if (iteration == 0)
  {
    reset_search (search);
    search->total_size = total;
    reallocate_search (search);

    if ((avoid >= 0) && (avoid < total))
      search->log[avoid] = TRUE;

    for (i = 0; i < seed_total; i++)
    {
      if ((seed[i] >= 0) && (seed[i] < total))
      {
        search->log[seed[i]] = TRUE;
        search->target[search->radius][search->total[search->radius]] = seed[i];
        search->origin[search->radius][search->total[search->radius]] = NEITHER;
        search->total[search->radius]++;
      }

      else
        exit (fprintf (global.outfile,
          "ERROR depth_search: seed value inappropriate (%d).\n", seed[i]));
    }
  }

  else
    search->radius--;

  if (search->radius >= 0)
  {
    for (terminus_found = FALSE; !terminus_found;)
    {
      node = search->target[search->radius][search->count[search->radius]];
      search->target[search->radius + 1][search->total[search->radius + 1]] =
        NEITHER;

      for
      (
        search->neighbor_id[search->radius] = 0;
        (neighbor = get_neighbor
          (nodes, node, search->neighbor_id[search->radius])) != EOF;
        search->neighbor_id[search->radius]++
      )
      {
        if (flag_neighbor)
          flag_neighbor
          (
            nodes,
            node,
            search->neighbor_id[search->radius],
            !search->log[neighbor]
          );

        if (!search->log[neighbor])
        {
          search->target
            [search->radius + 1]
            [search->total[search->radius + 1]] = neighbor;
          search->origin
            [search->radius + 1]
            [search->total[search->radius + 1]] = node;
          search->log[neighbor] = TRUE;

          break;
        }
      }

      if (search->target
        [search->radius + 1]
        [search->total[search->radius + 1]] == NEITHER)
        terminus_found = TRUE;

      else
      {
        search->count[search->radius + 1] = 
          search->total[search->radius + 1]++;
        search->radius++;
      }
    }

    return search->target[search->radius][search->count[search->radius]];
  }

  else
  {
    search->complete_flag = TRUE;
    return EOF;
  }
}


/* ///////////////////////////////////////////////////////////////// */

int get_search_origin
(
  SEARCH	*search,
  int		target,
  int		radius,
  int		count
)
{
  if ((target == NEITHER) && (radius == NEITHER) && (count == NEITHER))
    return search->origin[search->radius][search->count[search->radius]];

  else if (search->complete_flag)
  {
    if (target != NEITHER)
    {
      for (radius = 0; radius < search->total_size; radius++)
        for (count = 0; count < search->total[radius]; count++)
          if (search->target[radius][count] == target)
            return search->origin[radius][count];

      return NEITHER;
    }

    else if ((radius >= 0) && (radius < search->total_size) &&
      (count >= 0) && (count < search->total[radius]))
      return search->origin[radius][count];

    else
      return NEITHER;
  }

  else
    exit (fprintf (global.outfile,
      "ERROR get_search_target: search not complete.\n"));

  return TRUE;
}

/* ///////////////////////////////////////////////////////////////// */

int get_search_target
(
  SEARCH	*search,
  int		origin,
  int		radius,
  int		count
)
{
  if ((origin == NEITHER) && (radius == NEITHER) && (count == NEITHER))
    return search->target[search->radius][search->count[search->radius]];

  else if (search->complete_flag)
  {
    if ((origin >= 0) && (origin < search->total_size))
    {
      for (radius = 0; radius < search->total_size; radius++)
        for (count = 0; count < search->total[radius]; count++)
          if (search->origin[radius][count] == origin)
            return search->target[radius][count];

      return NEITHER;
    }

    else if ((radius >= 0) && (radius < search->total_size) &&
      (count >= 0) && (count < search->total[radius]))
      return search->target[radius][count];

    else if (count != NEITHER)
      return NEITHER;
  }

  else
    exit (fprintf (global.outfile,
      "ERROR get_search_target: search not complete.\n"));

  return TRUE;
}


/* ///////////////////////////////////////////////////////////////// */

int get_search_radius
(
  SEARCH	*search,
  int		target,
  int		origin
)
{
  int radius;
  int count;

  if ((origin == NEITHER) && (target == NEITHER))
    return search->radius;

  else if (search->complete_flag)
  {
    if ((target >= 0) && (target < search->total_size))
    {
      for (radius = 0; radius < search->total_size; radius++)
        for (count = 0; count < search->total[radius]; count++)
          if (search->target[radius][count] == target)
            return radius;

      return NEITHER;
    }

    if ((origin >= 0) && (origin < search->total_size))
    {
      for (radius = 0; radius < search->total_size; radius++)
        for (count = 0; count < search->total[radius]; count++)
          if (search->origin[radius][count] == origin)
            return radius;

      return NEITHER;
    }

    return NEITHER;
  }

  else
    exit (fprintf (global.outfile,
      "ERROR get_search_radius: search not complete.\n"));

  return TRUE;
}

