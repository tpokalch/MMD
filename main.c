
//PROBLEMS: define boltzman constant, MOL_MASS
//density = NUM_OF_MOL * MOL_MASS / BOX_DIMENSION ^ 3

#include <mlx.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define WIDTH 500
#define HEIGHT 500

#define NUM_OF_MOL 200
#define	BOX_DIMENSION WIDTH
#define MOL_SIZE 5 
#define MOL_MASS 1
#define k 1 
#define T 10 // Temperature in Kelvins
#define AVARAGE_SPEED sqrt(8 * k * T / (M_PI * MOL_MASS))
#define E0 50
#define DT 0.0001

float mymod(float x, float m)
{
        float ret;

        ret = fmod(x, m);
        if (x >= 0)
                return (ret);
        return (ret + m);
}

typedef struct	s_vector 
{
	float x;
	float y;
	float z;
}		t_vector;


typedef struct s_mol
{
	t_vector s;	// speed
	t_vector p;	// position
	t_vector a;	// acceleration
}		t_mol;

typedef struct		s_global
{
	void 			*mlx_ptr;
	void			*win_ptr;
	void			*img_ptr;
	int				*data_ptr;
	int				bpp;
	int				sz_l;
	int				e;


	t_mol			*mol;
}			t_global;


t_vector	get_random_unit_vector_in_box()
{
	t_vector ret;

	float theta = ( rand() / (float)RAND_MAX) * M_PI / 2.0;
	float phi = (rand() / (float)RAND_MAX) * M_PI / 2.0;

	ret.x = sin(theta) * cos(phi);	
	ret.y = sin(theta) * sin(phi);	
	ret.z = cos(theta);
	return (ret);
}

t_vector	get_random_unit_vector()
{
	t_vector ret;

	float theta = (rand() / (float)RAND_MAX) * (M_PI);
	float phi = (rand() / (float)RAND_MAX) * (M_PI * 2.0);
	ret.x = sin(theta) * cos(phi);	
	ret.y = sin(theta) * sin(phi);	
	ret.z = cos(theta);
	return (ret);

}

t_vector	scale(float s, t_vector v)
{
	t_vector ret;
	ret.x = v.x * s;
	ret.y = v.y * s;
	ret.z = v.z * s;
	return (ret);
}

t_vector	get_random_speed()
{
	int speed = mymod(rand(), (2 * AVARAGE_SPEED));
	t_vector unit_vector = get_random_unit_vector();
	t_vector ret = scale(speed, unit_vector);
/////////////
	ret.z = 0;
////////////
	return (ret);
}

t_vector diff(t_vector a, t_vector b)
{
	t_vector ret;

	ret.x = a.x - b.x;
	ret.y = a.y - b.y;
	ret.z = a.z - b.z;

	return (ret);
}

t_vector sum(t_vector a, t_vector b)
{
	t_vector ret;

	ret.x = a.x + b.x;
	ret.y = a.y + b.y;
	ret.z = a.z + b.z;

	return (ret);
}


float len(t_vector a)
{
	float ret = sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
	return (ret);
}

t_vector norm(t_vector a)
{
	float lent = len(a);
	a.x = a.x / lent;
	a.y = a.y / lent;
	a.z = a.z / lent;
	return (a);
}



float	dist(t_vector a, t_vector b)
{
	return (len(diff(a, b)));
}

void	init_random_pos(t_mol *mol, int i)
{
	printf("doing %d\n", i);

	int radius = rand();
	printf("radius is %d\n", radius);
	mol[i].p = get_random_unit_vector_in_box();
	printf("random unit vector: %f,%f,%f\n", mol[i].p.x, mol[i].p.y, mol[i].p.z);
	mol[i].p = scale(radius, mol[i].p);

	mol[i].p.x = mymod(mol[i].p.x, BOX_DIMENSION);
	mol[i].p.y = mymod(mol[i].p.y, BOX_DIMENSION);
	mol[i].p.z = 0/*mymod(mol[i].p.z, BOX_DIMENSION)*/;

	int j = 0;
	while (j < i)
	{
		if (dist(mol[j].p, mol[i].p) < 2 * MOL_SIZE)
		{
//			printf("mol %d has pos %f,%f,%f\n", j, mol[j].p.x, mol[j].p.y, mol[j].p.z); 
//			printf("mol %d has pos %f,%f,%f\n", i, mol[i].p.x, mol[i].p.y, mol[i].p.z);
			init_random_pos(mol, i);
		}
		j++;
	}
}

void	null_acceleration(t_mol *mol)
{
	int i = 0;
	while (i < NUM_OF_MOL)
	{
		mol[i].a.x = 0;
		mol[i].a.y = 0;
		mol[i].a.z = 0;
		i++;
	}
}

void	init_mol(t_mol *mol)
{
	srand(time(0));
	int i = 0;
	while (i < NUM_OF_MOL)
	{
		printf("initing mol %d\n", i);
//		check this
		(mol[i]).s = get_random_speed();
		init_random_pos(mol, i);
		printf("done %d\n", i);
//		null_acceleration(mol, i);
		i++;
	}
	null_acceleration(mol);
	printf(" %f percent done \n", (i + 1.0) / (float)NUM_OF_MOL);
}

float get_force(float r)
{
	float force = -48 * E0 / (float)MOL_SIZE * (pow(MOL_SIZE / r, 13) - 0.5 * pow(MOL_SIZE / r, 7));
//	if (r < 0 * MOL_SIZE)
//		force = 0;
	return(force);
}

void	write_acceleration(t_mol *mol, int i, int j, float distance)
{
	t_vector psti = mol[i].p;
	t_vector pstj = mol[j].p;

	t_vector from_j_to_i = diff(psti, pstj);
	from_j_to_i = norm(from_j_to_i);

	float force = get_force(distance);
//	printf("force on %d if %f\n", i, force);
	from_j_to_i = scale(force, from_j_to_i);
	t_vector from_i_to_j = scale(-1, from_j_to_i);
//	switched here
	mol[i].a = sum(mol[i].a, from_j_to_i);
	mol[j].a = sum(mol[j].a, from_i_to_j);
}

void	update_acceleration(t_mol *mol)
{
	float distance = 0;
	int i = NUM_OF_MOL - 1;
	while (i >= 0)
	{
		int j = 0;
		while (j < i)
		{
			if ((distance = dist(mol[i].p, mol[j].p)) < 2 * MOL_SIZE)
			{
				write_acceleration(mol, i, j, distance);
			}
			j++;
		}
		i--;
	}
}

void	update_speed(t_mol *mol)
{
	int i = 0;
	while (i < NUM_OF_MOL)
	{
		mol[i].s = sum(mol[i].s, scale(DT, mol[i].a));
		i++;
	}
}

void	update_pos(t_mol *mol)
{
	int i = 0;
	while (i < NUM_OF_MOL)
	{
		mol[i].p = sum(mol[i].p, scale(DT, mol[i].s));

		if (mol[i].p.x > BOX_DIMENSION || mol[i].p.x < 0)
			mol[i].p.x = mymod(mol[i].p.x, BOX_DIMENSION); 
		if (mol[i].p.y > BOX_DIMENSION || mol[i].p.y < 0)
			mol[i].p.y = mymod(mol[i].p.y, BOX_DIMENSION); 
		if (mol[i].p.z > BOX_DIMENSION || mol[i].p.z < 0)
			mol[i].p.z = mymod(mol[i].p.z, BOX_DIMENSION); 

		i++;
	}
}

float	dot(t_vector a, t_vector b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

void	draw_kinetic_energy(t_global *g, int x)
{
	float K = 0;
	int i = 0;

	while (i < NUM_OF_MOL)
	{
		K += dot(g->mol[i].s, g->mol[i].s);
		i++;
	}

	printf("drawing at x = %d\n", x);
	mlx_pixel_put(g->mlx_ptr, g->win_ptr, x, HEIGHT -lround(K), 0xFFFFFF);

//	mlx_pixel_put(g->mlx_ptr, g->win_ptr, x, lround(HEIGHT / 2.0), 0xFFFFFF);
	printf("K is %f\n", K);
}


int c = 0xFFFFFF;

void	draw_molecules(t_global *g)
{
//	c++;
	int i = 0;
	printf("c is %d\n", c);
	while (i < NUM_OF_MOL)
	{
		mlx_pixel_put(g->mlx_ptr, g->win_ptr, g->mol[i].p.x, HEIGHT - g->mol[i].p.y, c);
		i++;
	}

}

float	maxwell(float x)
{
	x = x / 50.0;
	float A = 20;
	float alpha = 0.1;
	float ret = A * x * x * exp(-alpha * x * x);
	return (ret);
}



void	draw_maxwell(t_global *g)
{
	static int flag = 0;
	flag++;

	int i = 0;

	while (i < WIDTH)
	{
		mlx_pixel_put(g->mlx_ptr, g->win_ptr, i, HEIGHT / 2 - maxwell(i), c);
		if (flag < 2)
		{
			if (i % 2)
			{
				mlx_pixel_put(g->mlx_ptr, g->win_ptr, i, HEIGHT / 2 - maxwell(i) -25 + (rand()) % 50, 0xFF0000);
			}
			if ((i + 1) % 6 == 1)
			{
				mlx_pixel_put(g->mlx_ptr, g->win_ptr, i, HEIGHT / 2 - maxwell(i) - 35 + (rand()) % 75, 0xFF0000);
			}
			else if ((i + 1) % 6 == 3)
			{
				mlx_pixel_put(g->mlx_ptr, g->win_ptr, i, HEIGHT / 2 - maxwell(i) - 55 + (rand()) % 100, 0xFF0000);
			}
			else if ((i + 1) % 6 == 5)
			{
				mlx_pixel_put(g->mlx_ptr, g->win_ptr, i, HEIGHT / 2 - maxwell(i) - 70 + (rand()) % 200, 0xFF0000);
			}


		}
		i++;
	}
}

int	loop(void *p)
{
        t_global *g = (t_global *)p;
	int i = 0;
	while (i < WIDTH)
	{
//		printf("position of %d is %f,%f,%f\n", 0, g->mol[0].p.x, g->mol[0].p.y, g->mol[0].p.z); 
		null_acceleration(g->mol);
		update_acceleration(g->mol);
		update_speed(g->mol);
		update_pos(g->mol);

//		mlx_clear_window(g->mlx_ptr, g->win_ptr);


//		draw_kinetic_energy(g, i);
		draw_molecules(g);
//		draw_maxwell(g);
		i++;
	}
	return (0);
}

int	main()
{
	t_global g;

	g.mlx_ptr = mlx_init();
	// with mymlx window should be created always after mlx_init
	// not necessery here, only for homogenity between versions
	g.win_ptr = mlx_new_window(g.mlx_ptr, WIDTH, HEIGHT, "window1");
	g.img_ptr = mlx_new_image(g.mlx_ptr, WIDTH, HEIGHT);
	g.data_ptr = (int *)mlx_get_data_addr(g.img_ptr, &g.bpp, &g.sz_l, &g.e);
//	g.win_ptr = mlx_new_window(g.mlx_ptr, WIDTH, HEIGHT, "window1");

	t_mol mol[NUM_OF_MOL];
	t_mol *a = (t_mol *)mol;
	g.mol = a;
	init_mol(a);



	mlx_loop_hook(g.mlx_ptr, loop, &g);
	mlx_loop(g.mlx_ptr);

//	loop(g.mlx_ptr);

	return (0);
}







